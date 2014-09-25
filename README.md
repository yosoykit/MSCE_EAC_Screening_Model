MSCE_EAC_Screening_Model
========================

These are the functions necessary to implement the MSCE-EAC Screening Method developed by Kit Curtius with contributions by Georg Luebeck.

`MSCE_EAC_Screening_Model.R`
- main source code, other functional scripts called from this file
- some user input for simulation size and age of patients are modifiable in this file
- Entire MSCE-EAC hybrid simulation algorithm and screening module performed through this file, data organized before screening implementation by script in `MSCE_EAC_pre_biop_setup.R`
- output provided for neoplastic prevalences based on biopsy-based or imaging based screening
- also computes probabilities of a missed malignancy in a dysplastic patient

`MSCE_EAC_screening_functions.R`
- contains all functions for cellular growth and computational screening implementation (biopsy-based or imaging)
- biopsy-based protocol may be adjusted to any protocol, forcep size, etc. Seattle protocol (standard) is currently provided
- Biopsy-based screening may occur through efficient and robust circular assumption on a plane. It may also occur through allowance of more diffusive clone growth (functions for which are in `MSCE_EAC_HexGrid.R`, `MSCE_EAC_neighborList.R`, `MSCE_EAC_genShape.R`)
- The geometric function for circular clones overlapping rectangular biopsies requires numerical integration, which we implement in this code via Legendre Gauss quadrature written in Fortran by Georg Luebeck and provided in this repo

`ex_parameter_list.R`
- example parameters provided for the biological rates and GERD prevalence values. Example given for all race males, fit to SEE incidence data and used for examples in Kong et al. 2014 and Screening paper.

`BE_density_func.R`
- function to compute rate of GERD-dependent BE onset rate
- outputs density and cumulative function for BE onset times assuming an exponential distribution (can be changed to any distrubtion desired)

`MSCE_EAC_incidence_projections.R`
- computes EAC incidence projections for sample BE cohort data produced by `MSCE_EAC_Screening_Model.R`
- performs ablation via fraction decimation controlled by matrix omega, which can be chosen by user based on desired efficacy of depleting cell types via radio-frequency ablation
