# SunflowerG
Sunflower rhizosphere microbiomes from 95 genotypes at Carrington, ND\
Analyses build on Pogoda et al. 2024 by including more community-level diversity analyses, analysis of multiple taxonomic levels, comparison of archaea/bacteria and fungi, and assessment of assembly mechanisms\
Analyses by Cliff Bueno de Mesquita, May 2024\
Finalized December 2024 and submitted to Phytobiomes Journal\
Revised May 2025\

## Repository Structure
In addition to the README, there are 4 subdirectories: code, data, FinalFigs, InitialFigs

code: contains 6 R scripts. The main analysis with all of the results for the paper is CarringtonGenotypeAnalysis_CB.R. The script is organized into different sections which can be viewed in RStudio outline. ReassignTaxonomy_CB.R is how taxonomy was reassigned in 2024 using the most recent SILVA and UNITE databases at the time. The other 4 scripts contain functions that are run from within the main CarringtonGenotypeAnalysis_CB.R script.

data: contains all of the input data needed to run the CarringtonGenotypeAnalysis_CB.R script

FinalFigs: folder containing the final figures for the manuscrip

InitialFigs: contains other exploratory figures that are not featured in the manuscript

