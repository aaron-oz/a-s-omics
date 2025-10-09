# purpose
this code runs a joint spatial transcriptomics model in inlabru to estimate:
 - the total feature counts with a poisson distribution
 - the feature density (counts per total counts) of a specific feature
 - the target feature counts (feature density * total counts)

it is meant to allow reproduction of the work and plots featured in:

*Spatial Modeling of Tissues for Morphogenic Field Analysis* as submitted for peer review oct 2025

# use
- download the repo
- set the repo location on the  01-array-run-job.R script
- source the 01-array-run-job.R script
- if/as needed, the 01 script could be launched in parrallel
