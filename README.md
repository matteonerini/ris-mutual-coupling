# Global Optimal Closed-Form Solutions for Intelligent Surfaces With Mutual Coupling: Is Mutual Coupling Detrimental or Beneficial?

This repository contains code to reproduce the paper:

> M. Nerini, H. Li, B. Clerckx, "[Global optimal closed-form solutions for intelligent surfaces with mutual coupling: Is mutual coupling detrimental or beneficial?](https://ieeexplore.ieee.org/document/11146418)," IEEE Trans. Wireless Commun., 2025.

## Code

The script `main_gen_ZIIs.m` generates the file `ZIIs-d-lambda-n.mat`, for *n* = 1, ..., 4, containing the mutual coupling matrices when the ratio between the wavelength and the inter-element distance is *n*.
These `.mat` files are then used by the other scripts to plot all the figures.

The script `main_Fign.m` reproduces Fig. *n*, for *n* = 3, ..., 6.
