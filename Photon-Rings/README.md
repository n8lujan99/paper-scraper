EHT Project Cheat Sheet
Title
Extracting Time Delays of Sgr A* from Light Curves

Authors / Colaborators
Dashon Jones, Prashant Kocherlakota, Paul Tiede, Dominic Chang, and Koushik Chatterjee, EHT

ABSTRACT
[What this project tries to answer, why it matters]
Time delay light curve analysis of photon ring images including multiple order emission. Using a plasma hotspot model orbiting
Sgr A*, light curves are generated with general relativistic ray-tracing to compare with EHT and ngEHT
observations. Using Gaussian processes and Bayesian inference, we are able to extract potential measurements
for the Lyapunov exponent for photon orbits. Specifically, we find the value of the Lyapunov to be π, consistent
with analytic predictions for Schwarzschild black holes, and the time delay values to match our theoretical
predictions.

Key words
[Concise set of topical tags]
General Relativity, SMBHs, Spin, Extreme Gravitational Lensing, Photon rings 
Time delay
very long baseline interferometry
Lyapunov exponent

Physics Topics to Review
[List core theory areas you need to brush up on]
⦁	Photon Rings, 
⦁	Gravitational Lensing
⦁	Zero Spin BH 
⦁	Null geodesics


Methods
[Key techniques, codes, datasets used]
THEMIS
Generate movies of plasma hotspots orbiting zero spin black holes and using the general relativistic, ray-traced, semi-analytic data to analyze light curves of various radii and observer inclination.
eht-imaging (Chael et al. (2018) --> VLBI data
uses Themis for interferometric visibilities

Results
[Main findings, important numbers, figures]


Next steps
[Where to push the project forward]




Open questions
[What’s unresolved, what you’re still skeptical about]
Assumes the black hole has zero spin
neglect shearing, expansion of hotspots scattering
optically thin



References
[Only essential ones to revisit quickly]
Broderick et al. 2022 - photon rings
Tiede, P., Pu, H.-Y., Broderick, A. E., et al. 2020, The Astrophysical Journal, 892, 132, doi: 10.3847/1538-4357/ab744c (Shearing and expansion)
Event Horizon Telescope Collaboration, Akiyama, K., Alberdi, A., et al. 2019, ApJL, 875, L1, doi: 10.3847/2041-8213/ab0ec7 —. 2022, ApJL, 930, L12, doi: 10.3847/2041-8213/ac6674 —. 2024, ApJL, 964, L26, doi: 10.3847/2041-8213/ad2df1
Walia, R. K., Kocherlakota, P., Chang, D. O., & Salehi, K. 2025, Physical Review D, 111, doi: 10.1103/physrevd.111.104074
Data: https://simbad.u-strasbg.fr/simbad/sim-ref?querymethod=bib&simbo=on&submit=submit+bibcode&bibcode=2019ApJ...875L...1E

Tasks
[Concrete to-dos related to this project]








Questions
[Things you want to ask advisor/collaborators, or figure out]




Helpful Resources
[Datasets, repos, lecture notes, textbooks, cheat links]

IPOLE User Guide for TACC 
https://docs.google.com/presentation/d/1ujYdDilKvcDeC1WKsLQQWuCmT1ux-lZC1WnQ321lUHg/edit?slide=id.p#slide=id.p

IPOLE User Guide Updated
https://docs.google.com/document/d/1-Vxn8ykR8qSwEljMDh_sAVkTiwIChJhaZa9hI_IXH-4/edit?tab=t.0

IPOLE Dump Files
https://app.box.com/s/udreucxy5kmru6mfcen3xv5ts70i6efr?page=1

Richard Computing Documentation 
https://docs.google.com/document/d/15DWakRSsxVBbCxMbJgBVTV839HIHuclXx9sT1OLTTcE/edit?tab=t.0

M87 Emission papers
https://drive.google.com/drive/u/0/folders/1URlu1oXRQ_4NbFBmBsp5waeNyWnyvkYS


On the Comparison of AGN with GRMHD Simulations: II. M87
Richard Anantua1,2,3,4,5 ⋆, Angelo Ricarte2,3, George Wong6,7, Razieh Emami,3 Roger Blandford4, Lani Oramas1, Hayley West1, Joaquin Duran1and Brandon Curd1,2,3
https://arxiv.org/pdf/2309.05602

The Pauli exclusion principle at strong coupling: Holographic matter and momentum space
Richard J. Anantua, Sean A. Hartnoll, Victoria L. Martin, David M. Ramirez
https://arxiv.org/pdf/1210.1590


Notes: 
Using time varrying pigeons julia script now. 
