using Plots 
using JuMP
using CPUTime
using Dates
using Random
using CPLEX
using Printf

# Définition de constantes pour le statut de résolution du problème
const OPTIMAL = JuMP.MathOptInterface.OPTIMAL
const INFEASIBLE = JuMP.MathOptInterface.INFEASIBLE
const UNBOUNDED = JuMP.MathOptInterface.DUAL_INFEASIBLE
const FEASIBLE_POINT = JuMP.MathOptInterface.FEASIBLE_POINT 


# Renvoie un tableau a 2 dimensions avec la distance entre les points
function distance(tabX, tabY)
    n = length(tabX)
    d = zeros(n,n)
    for i in 1:n
        for j in 1:n
            d[i, j] = sqrt((tabX[i] - tabX[j])^2 + (tabY[i] - tabY[j])^2)
        end
    end
    return d
end

# Fonction utilitaire pour la distance entre deux points
function dist_points(x1, y1, x2, y2)
    return sqrt((x1 - x2)^2 + (y1 - y2)^2)
end

# Fonction donné dans le sujet, pour lire le fichier de placement des points
function Lit_fichier_Placement(nom_fichier, tabX, tabY, f)
  n=0
  open(nom_fichier) do fic
        for (i,line) in enumerate(eachline(fic)) 
              lg = split(line," ")      
              if (i==1) # La première ligne est le nombre de villes 'n'
                  n= parse(Int,lg[1])
              else # Les lignes suivantes contiennent les données des villes
                 push!(tabX,parse(Float64,lg[3]))
                 push!(tabY,parse(Float64,lg[4]))   
                 push!(f,parse(Float64,lg[5]))               
              end              
        end    
  end
  return n
end

# Fonction donné dans le sujet, pour dessiner l'instance
function Dessine_Instance_Placement(nom_fichier, tabX, tabY)
    nom_fichier_pdf = replace(nom_fichier, ".flp" => "_instance.pdf")
    println("Création du fichier pdf de l'instance: ", nom_fichier_pdf)

    Plots.plot(tabX,tabY,seriestype = :scatter, legend = false, title="Répartition des villes") # Ajout d'un titre
    Plots.savefig(nom_fichier_pdf)
end

# Fonction donnée dans le sujet, pour dessiner la solution de placement des points
function Dessine_Solution_Placement(nom_fichier, n, tabX, tabY, S)

    X=Float64[]
    Y=Float64[] 
    
    # Génère un nom de fichier unique pour la solution (UFLP ou CFLP)
    nom_fichier_sol_pdf = replace(nom_fichier, ".flp" => "_sol.pdf")

    println("Création du fichier pdf de la solution: ", nom_fichier_sol_pdf)

    Plots.plot(tabX,tabY,seriestype = :scatter, legend = false)
    
    for i=1:n
       min_val=10e10
       minj_idx=0
       # Trouve la ville la plus proche de i dans S
       for j=1:n
          # Utilise 'dist_points' ici pour la correction de UndefVarError: dist
          if ( (S[j]==1) && (min_val > dist_points(tabX[i],tabY[i],tabX[j],tabY[j])) )
              min_val=dist_points(tabX[i],tabY[i],tabX[j],tabY[j])                          
              minj_idx=j             
          end
       end
       if (i!=minj_idx) # Évite de dessiner un segment d'une ville à elle-même si elle est son propre centre
          empty!(X)
          empty!(Y)         
          push!(X,tabX[i])
          push!(X,tabX[minj_idx])   
          push!(Y,tabY[i])
          push!(Y,tabY[minj_idx])             
          # dessine le segment reliant i à son centre
          Plots.plot!(X,Y, legend = false)
       end
    end
    
    Plots.savefig(nom_fichier_sol_pdf) # Ecrit la courbe créée à la ligne précédente dans un fichier .pdf
    
end

# FONCTIONS DE RÉSOLUTION DES PROBLÈMES D'OPTIMISATION

# Fonction pour résoudre le problème UFLP (Uncapacitated Facility Location Problem)
function UFLP(alpha, beta, tabX, tabY, f, n, nom_fichier_base)
    m = Model(CPLEX.Optimizer)
    tabDistance = distance(tabX, tabY)

    @variable(m, u[1:n], Bin) # Bin car u[i] = 1 si il y a un centre, 0 sinon
    @variable(m, v[1:n, 1:n], Bin)  # Bin car v[i][j] = 1 si le centre de la ville i dessert la ville j, 0 sinon

    @objective(m, Min,
        sum(f[j] * u[j] for j in 1:n) + # Coûts fixes d'ouverture des centres
        sum(alpha * tabDistance[i,j] * v[i,j] + beta * v[i,j] for i in 1:n, j in 1:n) # Coûts de transport et gestion
    )

    @constraint(m, c1[i in 1:n], sum(v[i,j] for j in 1:n) == 1) # Chaque ville est affectée à un centre
    @constraint(m, c2[i in 1:n, j in 1:n], v[i,j] <= u[j]) # Affectation seulement à un centre ouvert
    @constraint(m, c3[i in 1:n], v[i,i] == u[i]) # Si une ville est un centre, elle se dessert

    t1 = now()
    set_silent(m)   # pour cacher les logs
    optimize!(m)
    t2 = now()
    temps = t2 - t1

    status = termination_status(m)
    val = -1.0
    nb_centres = 0
    solution_u = nothing
    statut_str = "INFAISABLE"

    if status == OPTIMAL || status == FEASIBLE_POINT
        val = objective_value(m)
        solution_u = [value(u[i]) for i in 1:n]
        nb_centres = count(x -> x > 0.5, solution_u)
        
        if status == OPTIMAL
            statut_str = "OPTIMAL"
        else
            statut_str = "FAISABLE"
        end

        println("✅ Solution optimale trouvée : ")
        println("🏙️ Nombre total de villes : ", n)
        println("🏢 Nombre de centres ouverts : ", nb_centres)
        println("🎯 Objectif (coût total) : ", round(val, digits=2))
        println("⏱️ Temps de résolution : ", round(Millisecond(temps).value / 1000, digits=3), " secondes")

        Dessine_Solution_Placement(replace(nom_fichier_base, ".flp" => "_UFLP.flp"), n, tabX, tabY, solution_u)
    else
        println("❌ Pas de solution optimale. Statut : ", status)
    end

    return val, solution_u, nothing, statut_str
end

# Fonction pour résoudre le problème CFLP (Capacitated Facility Location Problem)
function CFLP(alpha, beta, tabX, tabY, f, n, q, nom_fichier_base)
    m = Model(CPLEX.Optimizer)
    tabDistance = distance(tabX, tabY)

    @variable(m, u[1:n], Bin)
    @variable(m, v[1:n, 1:n], Bin)

    @objective(m, Min,
        sum(f[j] * u[j] for j in 1:n) +
        sum(alpha * tabDistance[i,j] * v[i,j] + beta * v[i,j] for i in 1:n, j in 1:n)
    )

    @constraint(m, c1[i in 1:n], sum(v[i,j] for j in 1:n) == 1)                # Chaque ville est affectée
    @constraint(m, c2[i in 1:n, j in 1:n], v[i,j] <= u[j])                     # Affectation seulement à un centre ouvert
    @constraint(m, c3[i in 1:n], v[i,i] == u[i])                               # Si une ville est un centre, elle se dessert
    @constraint(m, c4[j in 1:n], sum(v[i,j] for i in 1:n) <= q * u[j])        # Contrainte de capacité

    t1 = now()
    set_silent(m)
    optimize!(m)
    t2 = now()
    temps = t2 - t1

    status = termination_status(m)
    val = -1.0
    nb_centres = 0
    solution_u = nothing
    statut_str = "INFAISABLE"

    if status == OPTIMAL || status == FEASIBLE_POINT
        val = objective_value(m)
        solution_u = [value(u[i]) for i in 1:n]
        nb_centres = count(x -> x > 0.5, solution_u)
        
        if status == OPTIMAL
            statut_str = "OPTIMAL"
        else
            statut_str = "FAISABLE"
        end

        println("✅ Solution optimale trouvée : ")
        println("🏙️ Nombre total de villes : ", n)
        println("🏢 Nombre de centres ouverts : ", nb_centres)
        println("🎯 Objectif (coût total) : ", round(val, digits=2))
        println("⏱️ Temps de résolution : ", round(Millisecond(temps).value / 1000, digits=3), " secondes")

        Dessine_Solution_Placement(replace(nom_fichier_base, ".flp" => "_CFLP.flp"), n, tabX, tabY, solution_u)
    else
        println("❌ Pas de solution optimale. Statut : ", status)
    end

    return val, solution_u, nothing, statut_str
end


# ALGORITHMES GLOUTONS

# Fonction gloutonne pour UFLP : affectation à base de distance uniquement
function glouton_uflp_distance(n, tabX, tabY, couts_fixes, alpha, nb_centres_a_ouvrir, val_optimale,nom_fichier_base)
    t1 = now()
    
    # Initialisation du statut par défaut
    statut_str = "INFAISABLE"
    val = -1.0
    solution_u = nothing
    solution_v = nothing
    
    if nb_centres_a_ouvrir > n
        println("❌ Erreur : Trop de centres demandés.")
        # Retourne directement avec INFAISABLE
    else
        # Calcul de la matrice des distances entre villes
        matrice_distance = distance(tabX, tabY)

        # Initialisation des variables de décision
        u_solution = zeros(Int, n)
        v_solution = zeros(Int, n, n)

        # Choix aléatoire des centres à ouvrir
        centres_ouverts = Random.randperm(n)[1:nb_centres_a_ouvrir]
        for j in centres_ouverts
            u_solution[j] = 1
        end

        # Pour chaque ville, on l'affecte au centre ouvert le plus proche
        for i = 1:n
            dist_min = Inf
            meilleur_centre = -1
            for j in centres_ouverts
                if matrice_distance[i, j] < dist_min
                    dist_min = matrice_distance[i, j]
                    meilleur_centre = j
                end
            end
            v_solution[i, meilleur_centre] = 1
        end

        # Calcul du coût total
        cout_total = 0.0
        for j = 1:n
            cout_total += couts_fixes[j] * u_solution[j]
        end
        for i = 1:n
            for j = 1:n
                cout_total += alpha * matrice_distance[i, j] * v_solution[i, j]
            end
        end

        # Vérification de la validité de la solution
        if cout_total != Inf && !isnan(cout_total) && cout_total >= 0
            val = cout_total
            solution_u = u_solution
            solution_v = v_solution
            nb_centres = sum(u_solution)
            
            # Détermination du statut basé sur la comparaison avec l'optimal
            if val_optimale > 0 && abs(val - val_optimale) < 1e-6
                statut_str = "OPTIMAL"
            else
                statut_str = "FAISABLE"
            end

            t2 = now()
            temps = t2 - t1
            
            println("✅ Solution glouton trouvée : ")
            println("🏙️ Nombre total de villes : ", n)
            println("🏢 Nombre de centres ouverts : ", nb_centres)
            println("🎯 Objectif (coût total) : ", round(val, digits=2))
            println("⏱️ Temps de résolution : ", round(Millisecond(temps).value / 1000, digits=3), " secondes")
            Dessine_Solution_Placement(replace(nom_fichier_base, ".flp" => "_GLOUTON_UFLP_distance.flp"), n, tabX, tabY, u_solution)
        else
            println("❌ Solution infaisable détectée. Coût : ", cout_total)
        end
    end

    return val, solution_u, solution_v, statut_str
end

# Heuristique gloutonne pour CFLP : affectation par proximité, en respectant les capacités
function glouton_cflp_distance(n, tabX, tabY, couts_fixes, alpha, nb_centres_a_ouvrir, capacite, val_optimale,nom_fichier_base)
    t1 = now()
    
    # Initialisation du statut par défaut
    statut_str = "INFAISABLE"
    val = -1.0
    solution_u = nothing
    solution_v = nothing
    
    if nb_centres_a_ouvrir > n
        println("❌ Erreur : Trop de centres demandés.")
        # Retourne directement avec INFAISABLE
    elseif nb_centres_a_ouvrir * capacite < n
        println("❌ Capacité totale insuffisante : ", nb_centres_a_ouvrir * capacite, " < ", n, " villes")
        # Retourne directement avec INFAISABLE
    else
        matrice_distance = distance(tabX, tabY)
        u_solution = zeros(Int, n)
        v_solution = zeros(Int, n, n)
        charge_centre = zeros(Int, n)
        solution_valide = true

        # Choix aléatoire des centres ouverts
        centres_ouverts = Random.randperm(n)[1:nb_centres_a_ouvrir]
        for j in centres_ouverts
            u_solution[j] = 1
            v_solution[j, j] = 1  # Le centre se sert lui-même
            charge_centre[j] += 1
            if charge_centre[j] > capacite
                println("❌ Dépassement de capacité initial")
                solution_valide = false
                break
            end
        end

        if solution_valide
            # Affectation des autres villes
            villes_non_affectees = [i for i in 1:n if u_solution[i] == 0]
            for i in villes_non_affectees
                dist_min = Inf
                meilleur_centre = -1
                for j in centres_ouverts
                    if charge_centre[j] < capacite && matrice_distance[i, j] < dist_min
                        dist_min = matrice_distance[i, j]
                        meilleur_centre = j
                    end
                end

                if meilleur_centre == -1
                    println("❌ Ville ", i, " non assignée : tous les centres sont pleins.")
                    solution_valide = false
                    break
                end

                v_solution[i, meilleur_centre] = 1
                charge_centre[meilleur_centre] += 1
            end
        end

        if solution_valide
            # Calcul du coût total
            cout_total = 0.0
            for j = 1:n
                cout_total += couts_fixes[j] * u_solution[j]
            end
            for i = 1:n
                for j = 1:n
                    cout_total += alpha * matrice_distance[i, j] * v_solution[i, j]
                end
            end

            # Vérification finale de la validité
            if cout_total != Inf && !isnan(cout_total) && cout_total >= 0
                val = cout_total
                solution_u = u_solution
                solution_v = v_solution
                nb_centres = sum(u_solution)
                
                # Détermination du statut basé sur la comparaison avec l'optimal
                if val_optimale > 0 && abs(val - val_optimale) < 1e-6
                    statut_str = "OPTIMAL"
                else
                    statut_str = "FAISABLE"
                end

                t2 = now()
                temps = t2 - t1
                
                println("✅ Solution glouton trouvée : ")
                println("🏙️ Nombre total de villes : ", n)
                println("🏢 Nombre de centres ouverts : ", nb_centres)
                println("🎯 Objectif (coût total) : ", round(val, digits=2))
                println("⏱️ Temps de résolution : ", round(Millisecond(temps).value / 1000, digits=3), " secondes")
                Dessine_Solution_Placement(replace(nom_fichier_base, ".flp" => "_GLOUTON_CFLP_distance.flp"), n, tabX, tabY, u_solution)
            else
                println("❌ Solution infaisable détectée. Coût : ", cout_total)
            end
        else
            println("❌ Impossible de construire une solution valide.")
        end
    end

    return val, solution_u, solution_v, statut_str
end


# ALGORITHMES GLOUTONS GÉOGRAPHIQUES

#FOnctions nécessaires pour les algorithmes gloutons géographiques
# Vérifie si une ville j est assignée
function est_assignee(j, v, n)
    somme = 0
    for i in 1:n
        somme += v[i, j]
    end
    return somme == 1
end

# Choisit la ville de départ (bien placée au centre)
function choisir_ville_depart(tabX, tabY, n)
    ville_choisie = 1
    meilleur_score = -Inf
    
    for i in 1:n
        # Score = position géographique + aléatoire pour diversité
        score = tabX[i] + tabY[i] + rand(-10:10)
        if score > meilleur_score
            meilleur_score = score
            ville_choisie = i
        end
    end
    
    return ville_choisie
end

# Compte le nombre de villes non assignées
function nb_villes_non_assignees(v, n)
    compteur = 0
    for j in 1:n
        if !est_assignee(j, v, n)
            compteur += 1
        end
    end
    return compteur
end

# Initialise la solution
function initialiser_solution(n)
    u = zeros(Int, n)
    v = zeros(Int, n, n)
    return u, v
end

# Calcule la distance vers le centre le plus proche
function distance_vers_centre_le_plus_proche(j, tabDistance, u, n)
    distance_min = Inf
    for centre in 1:n
        if u[centre] == 1
            if tabDistance[j, centre] < distance_min
                distance_min = tabDistance[j, centre]
            end
        end
    end
    return distance_min
end

# Choisit une ville éloignée avec de l'aléatoire
function choisir_ville_eloignee_avec_aleatoire(tabDistance, u, v, n)
    candidats_eloignes = Int[]
    distance_max = 0.0
    
    # Trouve la distance maximum
    for j in 1:n
        if !est_assignee(j, v, n)
            dist_min = distance_vers_centre_le_plus_proche(j, tabDistance, u, n)
            if dist_min > distance_max
                distance_max = dist_min
            end
        end
    end
    
    # Prend les candidats les plus éloignés (80% de la distance max)
    seuil = 0.8 * distance_max
    for j in 1:n
        if !est_assignee(j, v, n)
            dist_min = distance_vers_centre_le_plus_proche(j, tabDistance, u, n)
            if dist_min >= seuil
                push!(candidats_eloignes, j)
            end
        end
    end
    
    # Si aucun candidat, prendre une ville non assignée au hasard
    if isempty(candidats_eloignes)
        for j in 1:n
            if !est_assignee(j, v, n)
                push!(candidats_eloignes, j)
            end
        end
    end
    
    return rand(candidats_eloignes)
end

# Calcule la charge d'un centre
function calculer_charge_centre(centre, v, n)
    charge = 0
    for j in 1:n
        if v[centre, j] == 1
            charge += 1
        end
    end
    return charge
end

# Ouvre un centre
function ouvrir_centre!(ville_courante, u, v)
    u[ville_courante] = 1
    v[ville_courante, ville_courante] = 1  # Le centre se dessert lui-même
end

#Algorithmes gloutons géographiques (utilisent les fonctions précédentes)

# Algorithme glouton géographique pour UFLP
function glouton_geo_uflp(tabX, tabY, f, n, rayon_max, alpha, beta, nom_fichier_base, val_optimale)
    t1 = now()
    
    # Initialisation du statut par défaut
    statut_str = "INFAISABLE"
    val = -1.0
    solution_u = nothing
    solution_v = nothing
    
    # Calcul de la matrice des distances
    tabDistance = distance(tabX, tabY)

    # Initialisation
    u, v = initialiser_solution(n)

    # Choisir la première ville
    ville_courante = choisir_ville_depart(tabX, tabY, n)

    # Boucle principale avec protection contre boucle infinie
    max_iterations = n * 2
    iteration = 0
    
    while nb_villes_non_assignees(v, n) > 0 && iteration < max_iterations
        iteration += 1
        
        # Ouvrir un centre
        ouvrir_centre!(ville_courante, u, v)
        
        # Assigner les villes proches
        for j in 1:n
            if !est_assignee(j, v, n)
                if tabDistance[ville_courante, j] <= rayon_max
                    v[ville_courante, j] = 1
                end
            end
        end
        
        # Choisir la prochaine ville si nécessaire
        if nb_villes_non_assignees(v, n) > 0
            ville_courante = choisir_ville_eloignee_avec_aleatoire(tabDistance, u, v, n)
        end
    end

    # Vérification de la solution
    if nb_villes_non_assignees(v, n) == 0
        cout_objectif = calculer_objectif(u, v, tabDistance, f, n, alpha, beta)
        
        if cout_objectif != Inf && !isnan(cout_objectif) && cout_objectif >= 0
            val = cout_objectif
            solution_u = u
            solution_v = v
            nb_centres = sum(u)
            
            # Détermination du statut basé sur la comparaison avec l'optimal
            if val_optimale > 0 && abs(val - val_optimale) < 1e-6
                statut_str = "OPTIMAL"
            else
                statut_str = "FAISABLE"
            end

            t2 = now()
            temps = t2 - t1
            
            println("✅ Solution glouton géographique trouvée : ")
            println("🏙️ Nombre total de villes : ", n)
            println("🏢 Nombre de centres ouverts : ", nb_centres)
            println("🎯 Objectif (coût total) : ", round(val, digits=2))
            println("⏱️ Temps de résolution : ", round(Millisecond(temps).value / 1000, digits=3), " secondes")

            Dessine_Solution_Placement(replace(nom_fichier_base, ".flp" => "_GLOUTON_UFLP_GEO.flp"), n, tabX, tabY, u)
        else
            println("❌ Solution infaisable détectée. Objectif : ", cout_objectif)
        end
    else
        println("❌ Impossible d'assigner toutes les villes avec rayon_max = ", rayon_max)
        println("   Villes non assignées : ", nb_villes_non_assignees(v, n))
    end

    return val, solution_u, solution_v, statut_str
end

# Algorithme glouton géographique pour CFLP
function glouton_geo_cflp(tabX, tabY, f, n, q, rayon_max, alpha, beta, nom_fichier_base, val_optimale)
    t1 = now()
    
    statut_str = "INFAISABLE"
    val = -1.0
    solution_u = nothing
    solution_v = nothing

    # Calcul de la matrice des distances
    tabDistance = distance(tabX, tabY)
    
    # Initialisation
    u, v = initialiser_solution(n)
    
    # Choisir la première ville
    ville_courante = choisir_ville_depart(tabX, tabY, n)
    
    # Boucle principale
    while nb_villes_non_assignees(v, n) > 0
        # Ouvrir un centre
        ouvrir_centre!(ville_courante, u, v)
        
        # Assigner les villes proches avec limite de capacité
        charge_actuelle = calculer_charge_centre(ville_courante, v, n)
        capacite_restante = q - charge_actuelle
        nb_assignations = 0
        
        for j in 1:n
            if !est_assignee(j, v, n) && nb_assignations < capacite_restante
                if tabDistance[ville_courante, j] <= rayon_max
                    v[ville_courante, j] = 1
                    nb_assignations += 1
                end
            end
        end
        
        # Choisir la prochaine ville si nécessaire
        if nb_villes_non_assignees(v, n) > 0
            ville_courante = choisir_ville_eloignee_avec_aleatoire(tabDistance, u, v, n)
        end
    end
    
    if nb_villes_non_assignees(v, n) == 0
        cout_objectif = calculer_objectif(u, v, tabDistance, f, n, alpha, beta)
        
        if cout_objectif != Inf && !isnan(cout_objectif) && cout_objectif >= 0
            val = cout_objectif
            solution_u = u
            solution_v = v
            nb_centres = sum(u)
            
            # Détermination du statut basé sur la comparaison avec l'optimal
            if val_optimale > 0 && abs(val - val_optimale) < 1e-6
                statut_str = "OPTIMAL"
            else
                statut_str = "FAISABLE"
            end

            t2 = now()
            temps = t2 - t1
            
            println("✅ Solution glouton géographique trouvée : ")
            println("🏙️ Nombre total de villes : ", n)
            println("🏢 Nombre de centres ouverts : ", nb_centres)
            println("🎯 Objectif (coût total) : ", round(val, digits=2))
            println("⏱️ Temps de résolution : ", round(Millisecond(temps).value / 1000, digits=3), " secondes")
            
            Dessine_Solution_Placement(replace(nom_fichier_base, ".flp" => "_GLOUTON_CFLP_GEO.flp"), n, tabX, tabY, u)
        else
            println("❌ Solution infaisable détectée. Objectif : ", cout_objectif)
        end
    else
        println("❌ Impossible d'assigner toutes les villes avec les contraintes")
        println("   Villes non assignées : ", nb_villes_non_assignees(v, n))
        println("   Capacité par centre : ", q, ", Rayon max : ", rayon_max)
    end

    return val, u, v, statut_str
end

# Fonction pour calculer la valeur objective d'une solution
function calculer_objectif(u, v, tabDistance, f, n, alpha, beta)
    cout_fixe = sum(f[j] * u[j] for j in 1:n)
    cout_transport = sum(alpha * tabDistance[i,j] * v[i,j] + beta * v[i,j] for i in 1:n, j in 1:n)
    return cout_fixe + cout_transport
end


# META HEURISTIQUE

# Fonction pour générer une solution aléatoire avec des 0 et 1
function genere_solution_aleatoire(n, p)
    S = zeros(Int, n)
    indices_disponibles = collect(1:n)
    indices_choisis = []
    for i in 1:p
        idx_random = rand(1:length(indices_disponibles))
        push!(indices_choisis, indices_disponibles[idx_random])
        deleteat!(indices_disponibles, idx_random)
    end
    for i in indices_choisis
        S[i] = 1
    end
    return S
end

# Perturbe une solution en changeant la position d'un centre
function perturber_solution(n, S)
    S_new = copy(S)
    indices_avec_centre = findall(x -> x == 1, S_new)
    indices_sans_centre = findall(x -> x == 0, S_new)
    if !isempty(indices_avec_centre) && !isempty(indices_sans_centre)
        i = rand(indices_avec_centre)
        j = rand(indices_sans_centre)
        S_new[i], S_new[j] = S_new[j], S_new[i]
    end
    return S_new
end

# Calcule le coût total d'une solution
function cout_total_solution(n, tabX, tabY, f, S, alpha, beta, q = -1)
    tabDistance = distance(tabX, tabY)
    centres_ouverts = findall(x -> x == 1, S)
    if isempty(centres_ouverts)
        return Inf
    end
    v = zeros(Int, n, n)
    if q <= 0
        for i in 1:n
            couts = [alpha * tabDistance[i, j] + beta for j in centres_ouverts]
            meilleur_centre_idx = argmin(couts)
            meilleur_centre = centres_ouverts[meilleur_centre_idx]
            v[i, meilleur_centre] = 1
        end
    else
        charge_centres = zeros(Int, n)
        villes_couts = []
        for i in 1:n
            meilleur_cout = minimum([alpha * tabDistance[i, j] + beta for j in centres_ouverts])
            push!(villes_couts, (i, meilleur_cout))
        end
        sort!(villes_couts, by = x -> x[2])
        for (ville, _) in villes_couts
            meilleur_centre = -1
            meilleur_cout = Inf
            for centre in centres_ouverts
                if charge_centres[centre] < q
                    cout = alpha * tabDistance[ville, centre] + beta
                    if cout < meilleur_cout
                        meilleur_cout = cout
                        meilleur_centre = centre
                    end
                end
            end
            if meilleur_centre != -1
                v[ville, meilleur_centre] = 1
                charge_centres[meilleur_centre] += 1
            else
                return Inf
            end
        end
    end
    cout_fixe = sum(f[j] * S[j] for j in 1:n)
    cout_transport = sum(alpha * tabDistance[i,j] * v[i,j] + beta * v[i,j] for i in 1:n, j in 1:n)
    return cout_fixe + cout_transport
end

# Descente stochastique itérée pour UFLP
function descente_stochastique_iteree_uflp(tabX, tabY, f, n, p, alpha, beta, itermax, nbre_essais, nom_fichier_base, val_optimale)
    t1 = now()
    
    # Initialisation du statut par défaut
    statut_str = "INFAISABLE"
    val = -1.0
    solution_u = nothing
    solution_v = nothing
    
    meilleure_solution = genere_solution_aleatoire(n, p)
    meilleur_cout = cout_total_solution(n, tabX, tabY, f, meilleure_solution, alpha, beta, -1)
    
    for essai in 1:nbre_essais
        solution_actuelle = genere_solution_aleatoire(n, p)
        cout_actuel = cout_total_solution(n, tabX, tabY, f, solution_actuelle, alpha, beta, -1)
        for iter in 1:itermax
            solution_voisine = perturber_solution(n, solution_actuelle)
            cout_voisin = cout_total_solution(n, tabX, tabY, f, solution_voisine, alpha, beta, -1)
            if cout_voisin < cout_actuel
                solution_actuelle = solution_voisine
                cout_actuel = cout_voisin
            end
        end
        if cout_actuel < meilleur_cout
            meilleure_solution = solution_actuelle
            meilleur_cout = cout_actuel
        end
    end

    # Vérification de la solution
    if meilleur_cout != Inf && !isnan(meilleur_cout) && meilleur_cout >= 0
        val = meilleur_cout
        solution_u = meilleure_solution
        solution_v = nothing  # Pas besoin de v pour les métaheuristiques
        nb_centres = sum(meilleure_solution)
        
        # Détermination du statut basé sur la comparaison avec l'optimal
        if val_optimale > 0 && abs(val - val_optimale) < 1e-6
            statut_str = "OPTIMAL"
        else
            statut_str = "FAISABLE"
        end

        t2 = now()
        temps = t2 - t1
        
        println("✅ Solution descente stochastique trouvée : ")
        println("🏙️ Nombre total de villes : ", n)
        println("🏢 Nombre de centres ouverts : ", nb_centres)
        println("🎯 Objectif (coût total) : ", round(val, digits=2))
        println("⏱️ Temps de résolution : ", round(Millisecond(temps).value / 1000, digits=3), " secondes")
        
        Dessine_Solution_Placement(replace(nom_fichier_base, ".flp" => "_DESCENTE_UFLP.flp"), n, tabX, tabY, meilleure_solution)
    else
        println("❌ Aucune solution faisable trouvée. Meilleur coût : ", meilleur_cout)
    end
    
    return val, solution_u, solution_v, statut_str
end

# Descente stochastique itérée pour CFLP
function descente_stochastique_iteree_cflp(tabX, tabY, f, n, p, q, alpha, beta, itermax, nbre_essais, nom_fichier_base, val_optimale)
    t1 = now()
    
    # Initialisation du statut par défaut
    statut_str = "INFAISABLE"
    val = -1.0
    solution_u = nothing
    solution_v = nothing
    
    meilleure_solution = genere_solution_aleatoire(n, p)
    meilleur_cout = cout_total_solution(n, tabX, tabY, f, meilleure_solution, alpha, beta, q)
    
    # Augmentation automatique du nombre de centres si nécessaire
    while meilleur_cout == Inf && p < n
        p += 1
        meilleure_solution = genere_solution_aleatoire(n, p)
        meilleur_cout = cout_total_solution(n, tabX, tabY, f, meilleure_solution, alpha, beta, q)
    end

    if meilleur_cout != Inf
        for essai in 1:nbre_essais
            solution_actuelle = genere_solution_aleatoire(n, p)
            cout_actuel = cout_total_solution(n, tabX, tabY, f, solution_actuelle, alpha, beta, q)
            if cout_actuel == Inf
                continue
            end
            for iter in 1:itermax
                solution_voisine = perturber_solution(n, solution_actuelle)
                cout_voisin = cout_total_solution(n, tabX, tabY, f, solution_voisine, alpha, beta, q)
                if cout_voisin < cout_actuel
                    solution_actuelle = solution_voisine
                    cout_actuel = cout_voisin
                end
            end
            if cout_actuel < meilleur_cout
                meilleure_solution = solution_actuelle
                meilleur_cout = cout_actuel
            end
        end
    end

    # Vérification de la solution
    if meilleur_cout != Inf && !isnan(meilleur_cout) && meilleur_cout >= 0
        val = meilleur_cout
        solution_u = meilleure_solution
        solution_v = nothing
        nb_centres = sum(meilleure_solution)
        
        # Détermination du statut basé sur la comparaison avec l'optimal
        if val_optimale > 0 && abs(val - val_optimale) < 1e-6
            statut_str = "OPTIMAL"
        else
            statut_str = "FAISABLE"
        end

        t2 = now()
        temps = t2 - t1
        
        println("✅ Solution descente stochastique trouvée : ")
        println("🏙️ Nombre total de villes : ", n)
        println("🏢 Nombre de centres ouverts : ", nb_centres)
        println("🎯 Objectif (coût total) : ", round(val, digits=2))
        println("⏱️ Temps de résolution : ", round(Millisecond(temps).value / 1000, digits=3), " secondes")
        
        Dessine_Solution_Placement(replace(nom_fichier_base, ".flp" => "_DESCENTE_CFLP.flp"), n, tabX, tabY, meilleure_solution)
    else
        println("❌ Impossible de trouver une solution faisable même avec ", p, " centres")
    end
    
    return val, solution_u, solution_v, statut_str
end


#Fonctions pour les bornes inférieures relax, pour répondre à la question 7 
function borne_inferieure_relax_UFLP(alpha, beta, tabX, tabY, f, n)
    model = Model(CPLEX.Optimizer)
    tabDistance = distance(tabX, tabY)

    @variable(model, u[1:n] >= 0)
    @variable(model, v[1:n, 1:n] >= 0)

    @objective(model, Min,
        sum(f[j] * u[j] for j in 1:n) +
        sum(alpha * tabDistance[i,j] * v[i,j] + beta * v[i,j] for i in 1:n, j in 1:n)
    )

    @constraint(model, [i in 1:n], sum(v[i,j] for j in 1:n) == 1)
    @constraint(model, [i in 1:n, j in 1:n], v[i,j] <= u[j])
    @constraint(model, [i in 1:n], v[i,i] == u[i])

    set_silent(model)
    optimize!(model)

    if termination_status(model) == OPTIMAL
        return objective_value(model)
    else
        return Inf
    end
end

function borne_inferieure_relax_CFLP(alpha, beta, tabX, tabY, f, n, q)
    model = Model(CPLEX.Optimizer)
    tabDistance = distance(tabX, tabY)

    @variable(model, u[1:n] >= 0)
    @variable(model, v[1:n, 1:n] >= 0)

    @objective(model, Min,
        sum(f[j] * u[j] for j in 1:n) +
        sum(alpha * tabDistance[i,j] * v[i,j] + beta * v[i,j] for i in 1:n, j in 1:n)
    )

    @constraint(model, [i in 1:n], sum(v[i,j] for j in 1:n) == 1)
    @constraint(model, [i in 1:n, j in 1:n], v[i,j] <= u[j])
    @constraint(model, [i in 1:n], v[i,i] == u[i])
    @constraint(model, [j in 1:n], sum(v[i,j] for i in 1:n) <= q * u[j])

    set_silent(model)
    optimize!(model)

    if termination_status(model) == OPTIMAL
        return objective_value(model)
    else
        return Inf
    end
end


# MAIN
function main()
    nom_fichier_base = "inst_10000.flp"

    tabX = Float64[]
    tabY = Float64[]
    f = Float64[]

    # Lecture du fichier
    n = Lit_fichier_Placement(nom_fichier_base, tabX, tabY, f)
    println("Le fichier contient ",n, " villes")

    # Paramètres
    alpha = 1.0
    beta = 0.5
    q = 9 # Capacité maximale pour un centre dans le CFLP
    rayon_max = 50.0 # Rayon maximum pour l'algorithme glouton géographique
    nb_centres = max(1, min(n ÷ 15, Int(round(sqrt(n) * 0.8))))
    println("Paramètres : alpha = ", alpha, ", beta = ", beta, ", q = ", q, ", rayon_max = ", rayon_max, ", nb_centres = ", nb_centres)

    results = []

    # Résolutions exactes
    println("\n--- Résolution du problème UFLP exacte---")
    valUFLP, uUFLP, _, statusUFLP = UFLP(alpha, beta, tabX, tabY, f, n, nom_fichier_base)
    if uUFLP !== nothing
        push!(results, ("UFLP exact", sum(uUFLP), valUFLP, statusUFLP))
    else
        push!(results, ("UFLP exact", "---", valUFLP, statusUFLP))
    end
    
    println("\n--- Résolution du problème CFLP exacte---")
    valCFLP, uCFLP, _, statusCFLP = CFLP(alpha, beta, tabX, tabY, f, n, q, nom_fichier_base)
    if uCFLP !== nothing
        push!(results, ("CFLP exact", sum(uCFLP), valCFLP, statusCFLP))
    else
        push!(results, ("CFLP exact", "---", valCFLP, statusCFLP))
    end

    println("\n--- Résolution du problème UFLP (glouton) ---")
    valU_glouton, uU, _, statusU = glouton_uflp_distance(n, tabX, tabY, f, alpha, nb_centres, valUFLP,nom_fichier_base)
    push!(results, ("UFLP glouton", uU !== nothing ? sum(uU) : "---", valU_glouton, statusU))

    println("\n--- Résolution du problème CFLP (glouton) ---")
    nb_centres_CFLP = nb_centres * q < n ? ceil(Int, n / q) : nb_centres
    valC_glouton, uC, _, statusC = glouton_cflp_distance(n, tabX, tabY, f, alpha, nb_centres_CFLP, q, valCFLP,nom_fichier_base)
    push!(results, ("CFLP glouton", uC !== nothing ? sum(uC) : "---", valC_glouton, statusC))

    println("\n--- Résolution du problème UFLP (glouton géographique) ---")
    val_glouton_uflp, u_glouton_geo_uflp, _, status_glouton_geo_uflp = glouton_geo_uflp(tabX, tabY, f, n, rayon_max, alpha, beta, nom_fichier_base, valUFLP)
    push!(results, ("UFLP glouton géo", u_glouton_geo_uflp !== nothing ? sum(u_glouton_geo_uflp) : "---", val_glouton_uflp, status_glouton_geo_uflp))

    println("\n--- Résolution du problème CFLP (glouton géographique) ---")
    val_glouton_cflp, u_glouton_geo_cflp, _, status_glouton_geo_cflp = glouton_geo_cflp(tabX, tabY, f, n, q, rayon_max, alpha, beta, nom_fichier_base, valCFLP)
    push!(results, ("CFLP glouton géo", u_glouton_geo_cflp !== nothing ? sum(u_glouton_geo_cflp) : "---", val_glouton_cflp, status_glouton_geo_cflp))

    println("\n--- Résolution du problème UFLP (descente stochastique itérée) ---")
    p_estime = max(1, min(n ÷ 15, Int(round(sqrt(n) * 0.8))))
    val_descente_uflp, u_descente_uflp, _, status_descente_uflp = descente_stochastique_iteree_uflp(tabX, tabY, f, n, p_estime, alpha, beta, 1000, 75, nom_fichier_base, valUFLP)
    push!(results, ("UFLP descente", u_descente_uflp !== nothing ? sum(u_descente_uflp) : "---", val_descente_uflp, status_descente_uflp))

    println("\n--- Résolution du problème CFLP (descente stochastique itérée) ---")
    val_descente_cflp, u_descente_cflp, _, status_descente_cflp = descente_stochastique_iteree_cflp(tabX, tabY, f, n, p_estime, q, alpha, beta, 1000, 75, nom_fichier_base, valCFLP)
    push!(results, ("CFLP descente", u_descente_cflp !== nothing ? sum(u_descente_cflp) : "---", val_descente_cflp, status_descente_cflp))

    println("\n📊 COMPARAISON COMPLÈTE DES RÉSULTATS")
    println("Méthode           | Centres ouverts | Coût total   | Statut")
    println("------------------|-----------------|--------------|-----------")
    for (methode, nb, cout, statut) in results
        nb_str = string(nb)
        cout_str = (cout == -1.0 || cout == Inf || cout == nothing) ? "---" : @sprintf("%.2f", cout)
        println(rpad(methode,18), " | ", lpad(nb_str,15), " | ", lpad(cout_str,12), " | ", statut)
    end

    println("\n===========================================")
    println("📈 GAPS PAR RAPPORT À LA RELAXATION LINEAIRE")
    println("===========================================")
    println("gap = (zbestsol - zb) / zbestsol")

    println("Le gap mesure l'écart relatif entre une solution heuristique et une borne inférieure (relaxation linéaire).")
    println("Plus le gap est faible, plus la solution heuristique est proche de l'optimal théorique.")
    println()


    zb_UFLP = borne_inferieure_relax_UFLP(alpha, beta, tabX, tabY, f, n)
    zb_CFLP = borne_inferieure_relax_CFLP(alpha, beta, tabX, tabY, f, n, q)

    if statusUFLP == "OPTIMAL" || statusUFLP == "FAISABLE"
        println("\n🔍 Analyse UFLP :")
        println("   Borne inférieure UFLP (relax) : ", round(zb_UFLP, digits=2))

        if valU_glouton > 0
            gap = ((valU_glouton - zb_UFLP) / valU_glouton) * 100
            println("   📊 Gap UFLP (glouton) : ", round(gap, digits=2), "%")
            evaluer_performance(gap, "UFLP glouton", "la borne relaxée")
        end

        if val_glouton_uflp > 0
            gap = ((val_glouton_uflp - zb_UFLP) / val_glouton_uflp) * 100
            println("   📊 Gap UFLP (glouton géo) : ", round(gap, digits=2), "%")
            evaluer_performance(gap, "UFLP glouton géo", "la borne relaxée")
        end

        if val_descente_uflp != Inf
            gap = ((val_descente_uflp - zb_UFLP) / val_descente_uflp) * 100
            println("   📊 Gap UFLP (descente) : ", round(gap, digits=2), "%")
            evaluer_performance(gap, "UFLP descente", "la borne relaxée")
        end
    end

    if statusCFLP == "OPTIMAL" || statusCFLP == "FAISABLE"
        println("\n🔍 Analyse CFLP :")
        println("   Borne inférieure CFLP (relax) : ", round(zb_CFLP, digits=2))

        if valC_glouton > 0
            gap = ((valC_glouton - zb_CFLP) / valC_glouton) * 100
            println("   📊 Gap CFLP (glouton) : ", round(gap, digits=2), "%")
            evaluer_performance(gap, "CFLP glouton", "la borne relaxée")
        end

        if val_glouton_cflp > 0
            gap = ((val_glouton_cflp - zb_CFLP) / val_glouton_cflp) * 100
            println("   📊 Gap CFLP (glouton géo) : ", round(gap, digits=2), "%")
            evaluer_performance(gap, "CFLP glouton géo", "la borne relaxée")
        end

        if val_descente_cflp != Inf
            gap = ((val_descente_cflp - zb_CFLP) / val_descente_cflp) * 100
            println("   📊 Gap CFLP (descente) : ", round(gap, digits=2), "%")
            evaluer_performance(gap, "CFLP descente", "la borne relaxée")
        end
    end

    Dessine_Instance_Placement(nom_fichier_base, tabX, tabY)
end

function evaluer_performance(gap, nom_methode, contexte="")
    if gap <= 0.1
        println("      ✅✅ Parfait : ", nom_methode, " trouve la même valeur que l'optimal !!!")
    elseif gap < 3
        println("      ✅ Excellent : ", nom_methode, " est très proche de ", contexte)
    elseif gap < 7
        println("      👍 Bien : ", nom_methode, " donne une bonne solution avec ", contexte)
    elseif gap < 15
        println("      🔶 Moyen : ", nom_methode, " reste acceptable avec ", contexte)
    else
        println("      ⚠️ Limité : ", nom_methode, " peine avec ", contexte)
    end
end

#Appel de la fonction main pour démarrer le programme
main()