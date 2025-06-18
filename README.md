# Facility Location Optimization – UFLP / CFLP  
# Optimisation de placement de centres logistiques – UFLP / CFLP

## 🇬🇧 English version

This project was developed as part of a university module on combinatorial optimization. It focuses on modeling and solving two classic facility location problems:

- **UFLP** – Uncapacitated Facility Location Problem  
- **CFLP** – Capacitated Facility Location Problem

We implemented and compared several approaches:
- Exact resolution via MILP (JuMP + CPLEX)
- Greedy heuristics (distance-based rules)
- Stochastic local descent metaheuristic

### Objectives
- Model real-world logistics problems
- Compare exact vs approximate methods (in time and quality)
- Evaluate performance on different instance sizes

### Tools & Libraries
- Language: Julia
- Modeling: JuMP
- Solver: CPLEX
- Visualization: PyPlot / Matplotlib

### Repository Structure
- `projet.jl` : source code 
- `instances/` : Test datasets
- `RapportOptimisation.pdf` : Final project report (in French)

---

## 🇫🇷 Version française

Ce projet a été réalisé dans le cadre d’un module d’optimisation combinatoire. Il porte sur la modélisation et la résolution de deux problèmes classiques de placement de centres logistiques :

- **UFLP** – Problème de localisation sans contrainte de capacité  
- **CFLP** – Problème de localisation avec contrainte de capacité

Nous avons comparé différentes approches :
- Résolution exacte via PLNE (JuMP + CPLEX)
- Algorithmes gloutons (basés sur la distance)
- Métaheuristique de descente locale stochastique

### Objectifs
- Modéliser un problème logistique réaliste
- Comparer méthodes exactes et approchées (temps et qualité)
- Tester la montée en charge sur des instances variées

### Technologies utilisées
- Langage : Julia
- Modélisation : JuMP
- Résolution : CPLEX
- Visualisation : PyPlot / Matplotlib


### Contenu du dépôt
- `projet.jl` : Code source
- `RapportOptimisation.pdf` : Rapport final du projet (en français)

---

## 👤 Authors / Auteurs
- Catherine CALI
- Katell GOUZERH
