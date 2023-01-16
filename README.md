# SocialSecurityDesign
Code repository for the paper "Social Security Design and its Political Support". Programs are written in Matlab (R2018b) and runs as is (no specific dependencies). The Live Script ```EGM.mlx``` contains the main steps for producing quantitative results. 

## Changing the calibration of the Argentina case:

The code for the Argentina case runs as follows:
1. The file ```Inputdata.xlsx``` contains data on $\nu_t$ in the worksheet "PopulationData". 
2. The file ```EGM.mlx``` runs the main simulations.
  2.1. The targets for the calibration are defined in the settings section - in the struct ```target```:
  <img src="snippets/settings.png" height = "128">
