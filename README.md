# Metabolic modelling of Commensalibacter sp.

Sequence a Genome II

Authors: Perrine Steffe, Karim Saied, Boris Schnider

Supervised by Noushin Hadadi

The project aims at building a Flux Balance Model (FBM) representing a simplification of the Commensalibacter sp. metabolism, from the previously identified genes and reactions in this organism. The reconstructed model will allow to investigate, through growth simulations with different input metabolites, the minimal glucose and oxygen uptake rates required for a Commensalibacter sp. cell to grow.

This directory contains the files (script and models) required to run the Flux Balance Analysis for Commensalibacter sp .

To perform the analysis, make sure you have all the Dependencies installed on a UNIX system and follow the steps in the Pipeline section.

### Dependencies

- MATLAB R2017b (9.3.0.713579) with the following packages installed:
    - The Cobra Toolbox
    - The Raven Toolbox

### Script

- biomass_model.m : script that defines the biomass reaction with the Biomass Building Blocks (BBBs) found in Gluconobacter diazotrophicus model. It also contains all the code to uptake the metabolites required and fill the gaps, i.e add the missing reactions for biomass production, enabling the production of BBBs.

### Data

- glc_dia.mat : Gluconobacter diazotrophicus mis-annotated model from BioModels Database EMBL-EBI (“ BMID000000141568.xml ”)

- model_commensalibacter.mat : First Commensalibacter sp . model built from scratch with “The Raven toolbox” by Hadadi Noushin

- model21.mat : Final FBM, for which the biomass equation is solved to visualize the growth prediction

### Pipeline

1. Download the packages “The Cobra Toolbox” and “The Raven Toolbox” and drag them in your MATLAB folder

2. Add “cobratoolbox” and “RAVEN_functions” folders to your path in MATLAB

3. Load glc_dia.m and model_commensalibacter.mat in your MATLAB environment (by double-clicking on your models from MATLAB) or by using the command load( 'Model_name' );

4. Run the script biomass_model.m

5. The presence or absence of growth, as well as its extent, can be observed in sol_FBA.obj

6. Running the following lines with different uptaken compounds KEGG IDs (such as C00031 and C00007 for glucose and oxygen respectively) and values after the “=” sign allows to test the effect of changing the uptaken rates of specific metabolites on the growth of Commensalibacter sp. : model21.lb(find(ismember(model21.rxns, 'EXC_BOTH_C00031' )))=0; sol_FBA=optimizeCbModel(model21);
