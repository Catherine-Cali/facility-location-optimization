# Optimisation de placement de centres logistiques – UFLP / CFLP

Ce projet a été réalisé dans le cadre d’un module d’optimisation combinatoire. Il vise à modéliser et résoudre deux variantes classiques de problème de localisation de centres logistiques :

- **UFLP** (Uncapacitated Facility Location Problem)  
- **CFLP** (Capacitated Facility Location Problem)

Nous y avons comparé différentes méthodes :
- Résolution exacte via PLNE (avec JuMP + CPLEX)
- Algorithmes gloutons (basés sur la distance)
- Métaheuristique de descente locale stochastique

## Objectifs
- Comprendre les enjeux opérationnels de la localisation de centres
- Implémenter et comparer des méthodes exactes et approchées
- Évaluer la qualité et la rapidité des algorithmes selon la taille des instances

## Technologies utilisées
- **Langage :** Julia
- **Modélisation :** JuMP
- **Résolution :** CPLEX
- **Visualisation :** PyPlot / Matplotlib

## Pourquoi Julia ?
Ce projet a aussi été l’occasion de découvrir Julia, un langage performant pour l’optimisation, tout en bénéficiant de fonctions de chargement de données fournies par l’enseignant. Cela a simplifié le développement et accéléré les expérimentations.

## Contenu du dépôt
- `src/` : Code source des méthodes exactes et heuristiques
- `instances/` : Jeux d’instances testés
- `résultats/` : Comparaison des performances
- `rapport.pdf` : Rapport final du projet

## Auteurs
- Catherine CALI
- Katell GOUZERH
