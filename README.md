# Global Optimal Closed-Form Solutions for Intelligent Surfaces With Mutual Coupling: Is Mutual Coupling Detrimental or Beneficial?

This code package is related to the paper:

M. Nerini, H. Li, and B. Clerckx, "[Global Optimal Closed-Form Solutions for Intelligent Surfaces With Mutual Coupling: Is Mutual Coupling Detrimental or Beneficial?](https://ieeexplore.ieee.org/document/11146418)," IEEE Trans. Wireless Commun., 2025.

## Content of Code Package

The file `main_gen_ZIIs.m` generates `ZIIs-d-lambda-i.mat`, for *i* = 1, ..., 4, containing the mutual coupling matrices when the ratio between the wavelength and the inter-element distance is *i*.
These `.mat` files are then used by the other scripts to plot all figures.

The file `main_Figi.m` reproduces Fig. *i* in the paper, for *i* = 3, ..., 6.
