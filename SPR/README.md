# README

This folder contains code for calculating exact distances in ranked SPR space (HSPR and RSPR) on small trees (up to 7 leaves, rankedspr_distances.py).
We further implement ideas for approximating distances (rankedspr_approximations.py), and methods for testing them (rankedspr_testing_distances.py).
In rankedspr_exploation.py we analyse characteristics of ranked SPR spaces, including sizes of orbits, diameters, etc.

For calculating exact distances in RSPR and HSPR, we use the implementation of Seidel's algorithm (seidel.c) of @kelmes.