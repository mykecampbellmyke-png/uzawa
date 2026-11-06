# uzawa

> Implementation en R de l'algorithme d'Uzawa pour les problemes
> d'optimisation sous contraintes lineaires et non lineaires.

## Description

L'algorithme d'Uzawa est une methode iterative pour resoudre des
problemes de point selle, qui apparaissent naturellement en
optimisation sous contraintes et dans la resolution numerique
d'equations aux derivees partielles (par exemple le systeme de
Stokes en mecanique des fluides).

Ce package fournit une implementation generale qui fonctionne pour :

- une fonction objectif **lineaire ou non lineaire** ;
- des contraintes d'**egalite** et/ou d'**inegalite** ;
- des contraintes **lineaires ou non lineaires**.

Le principe est simple : a chaque iteration, on alterne

1. une **etape primale** : minimisation du Lagrangien sans contrainte ;
2. une **etape duale** : mise a jour des multiplicateurs, avec
   projection sur `R_+` pour les contraintes d'inegalite.

## Installation

Le package n'est pas sur le CRAN. On l'installe depuis GitHub :

```r
# install.packages("remotes")
remotes::install_github("mykecampbellmyke-png/uzawa")
