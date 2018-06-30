glc_dia = readCbModel('BMID000000141568.xml','fileType','SBML')


printRxnFormula(glc_dia)
glc_dia.mets=glc_dia.metNames; %    Replacementof the metabolites (mets) by 
%   the met names in the reactions
%   We need to rerun the printRxnFormula


%mets=model_commensalibacter.mets;
%rxns=model_commensalibacter.rxns;
%v_cat = vertcat(mets, rxns);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Extraction of the biomass building blocks (BBBs) in the Gluconobacter model
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

BBBs=strfind(glc_dia.mets,'bm') %   Find index in glc_dia containing extension “bm” 
index = cellfun('isempty',BBBs); %  Remove the index 0 in BBBs
biomass_names=glc_dia.metNames(~index); %   Extraction of the names corresponding to index of BBBs in glc_dia


%   Find the biomass name in the Commensalibacter model metNames
[a,b]=ismember(model_commensalibacter.metNames,biomass_names);


%   Manually, the mets names from biomass_names are converted into their corresponding
%   KEGG IDs and copy-pasted into biomass_kegg.
biomass_kegg = {};
[e,f] = ismember(model_commensalibacter.mets,biomass_kegg);
index_commens = find(f); %   Index of biomass_kegg in the Commensalibacter model are stored in index_commens.

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Addition of the biomasse reaction in the Commensalibacter model 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   The instance of the biomass reaction is created through the addition of
%   a met called Biomass_c.
metsToAdd.mets={'Biomass_c'};
metsToAdd.metNames={'Biomass'};
metsToAdd.metFormulas={'NA'};
metsToAdd.b=[0];
metsToAdd.compartments={'s'};

model_biomass=addMets(model_commensalibacter,metsToAdd);

%   Addition of the biomass reaction in the model
rxnsToAdd.rxns={'Biomass_c'};
rxnsToAdd.equations={'34.7965 C00001 + 40.1701 C00002 + 0.005623 C00003 + 0.005623 C00006 + 0.005623 C00010 + 0.005623 C00016 + 0.005623 C00019 + 0.25601 C00025 + 0.5958 C00037 + 0.50006 C00041 + 0.2091 C00044 + 0.33355 C00047 + 0.23468 C00049 + 0.28828 C00062 + 0.12988 C00063 + 0.25601 C00064 + 0.2097 C00065 + 0.005623 C00068 + 0.14934 C00073 + 0.14027 C00075 + 0.055157 C00078 + 0.18056 C00079 + 0.13425 C00082 + 0.08898 C00097 + 0.005623 C00101 + 0.43866 C00123 + 0.023503 C00131 + 0.092623 C00135 + 0.21543 C00148 + 0.23468 C00152 + 0.4116 C00183 + 0.24665 C00188 + 0.005623 C00234 + 0.005623 C00255 + 0.28255 C00407 + 0.023503 C00459 + 0.092476 C04574 + 0.005623 C00229 <=> 40 C00008 + 39.9944 C00009 + 0.60238 C00013 + 40 C00080 + 1 Biomass_c'};
rxnsToAdd.rxnNames={'Biomass reaction'};
rxnsToAdd.eccodes={'NA'};
rxnsToAdd.rxnComps=[1];
rxnsToAdd.lb=[0];
rxnsToAdd.ub=[50];
rxnsToAdd.c=[1];

model4=addRxns(model_biomass,rxnsToAdd,1);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Definition of the extracellular environment
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   The solver is set %To tell cobra which solver(procedure) to use to solve 
%   the following line. LP = linear programming
changeCobraSolver('matlab','LP');  
sol_FBA=optimizeCbModel(model4); %  A first FBA is run to check whether biomass is produced.


%   Addition of common mets normally uptaken by bacteria
mets_to_uptake={};
[model5 addedRxns]=addExchangeRxns(model4,'BOTH',mets_to_uptake);


model5.lb(model5.lb==-1000)=-50; %  Only instance with -1000 are changed.
model5.ub(model5.ub==1000)=50; %  Only instance with 1000 are changed.


sol_FBA=optimizeCbModel(model5); %   The model is rerun

[BBB, BBBprod]=checkBBBTasks(model5); % The production of BBBs are checked.


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Curation of the model (gap filling, additional uptakes and secretion)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   uptake to produce L-Cysteine and dTTP
sulfide_to_uptake={'C00283'}; %   Sulfide uptake
[model6 addedRxns]=addExchangeRxns(model5, 'BOTH', sulfide_to_uptake);

%   Adding reaction to produce NAD+ and NADP+
rxnsToAdd.rxns={'R00481'}; %    Entry ID kegg pathways
rxnsToAdd.equations={'C00049 + C00007 <=> C05840 + C00027'} % Reaction formula
rxnsToAdd.rxnNames={'L-aspartate:oxygen oxidoreductase'}; % Reaction name
rxnsToAdd.eccodes={'1.4.3.16'}; % Reaction ID EC codes
rxnsToAdd.rxnComps=[1];
rxnsToAdd.lb=[-50]; %   Quantity of reagent consumed by the reaction
rxnsToAdd.ub=[50]; %    Quantity of products produced by the reaction
rxnsToAdd.c=[0];
 
model7=addRxns(model6,rxnsToAdd,1); %   Add reaction to the previous model.
 
[BBB, BBBprod]=checkBBBTasks(model7); %   Control of BBBs production and all metabolists according to the parameter
sol_FBA=optimizeCbModel(model7); %   Solve FBA of the new model to see if the bacteria is able to grow
 
%   Adding reaction to produce CoA
rxnsToAdd.rxns={'R00489'};
rxnsToAdd.equations={'C00049 <=> C00099 + C00011'}
rxnsToAdd.rxnNames={'L-aspartate 1-carboxy-lyase (beta-alanine-forming)'};
rxnsToAdd.eccodes={'4.1.1.11'};
rxnsToAdd.rxnComps=[1];
rxnsToAdd.lb=[-50];
rxnsToAdd.ub=[50];
rxnsToAdd.c=[0];
 
model8=addRxns(model7,rxnsToAdd,1);
 
rxnsToAdd.rxns={'R02472'};
rxnsToAdd.equations={'C00522 + C00006 <=> C00966 + C00005 + C00080'}
rxnsToAdd.rxnNames={'(R)-Pantoate:NADP+ 2-oxidoreductase'};
rxnsToAdd.eccodes={'1.1.1.169'};
rxnsToAdd.rxnComps=[1];
rxnsToAdd.lb=[-50];
rxnsToAdd.ub=[50];
rxnsToAdd.c=[0];
 
model9=addRxns(model8,rxnsToAdd,1);
[BBB, BBBprod]=checkBBBTasks(model9);
 
%   Adding a reaction in the pyrimidine pathway to produce CTP, UDP and TTP
rxnsToAdd.rxns={'R01869'};
rxnsToAdd.equations={'C00337 + C00003 <=> C00295 + C00080 + C00004'}
rxnsToAdd.rxnNames={'(S)-dihydroorotate:NAD+ oxidoreductase'};
rxnsToAdd.eccodes={'1.3.1.14'};
rxnsToAdd.rxnComps=[1];
rxnsToAdd.lb=[-50];
rxnsToAdd.ub=[50];
rxnsToAdd.c=[0];

model10=addRxns(model9,rxnsToAdd,1);
[BBB, BBBprod]=checkBBBTasks(model10);
printRxnFormula(model10);
 
%   Adding reaction to produce TTP
rxnsToAdd.rxns={'R02094'};
rxnsToAdd.equations={'C00002 + C00364 <=> C00008 + C00363'}
rxnsToAdd.rxnNames={'ATP:dTMP phosphotransferase'};
rxnsToAdd.eccodes={'2.7.4.9'};
rxnsToAdd.rxnComps=[1];
rxnsToAdd.lb=[-50];
rxnsToAdd.ub=[50];
rxnsToAdd.c=[0];
 
model11=addRxns(model10,rxnsToAdd,1);
[BBB, BBBprod]=checkBBBTasks(model11);
 
%   Adding reaction to produce Tetrahydrofolate and 10-Formyltetrahydrofolate
metsToAdd.mets={'C21615'}; %   Glycolaldehyde triphosphate
metsToAdd.metNames={'Glycolaldehyde_triphosphate'};
metsToAdd.metFormulas={'C2H7O11P3'};
metsToAdd.b=[0];
metsToAdd.compartments={'s'};
model12 = addMets(model11,metsToAdd);
 
rxnsToAdd.rxns={'R11719'};
rxnsToAdd.equations={'C04895 <=> C01300 + C21615'}
rxnsToAdd.rxnNames={'7,8-dihydroneopterin-3-triphosphate glycolaldehyde-phosphate-lyase'};
rxnsToAdd.eccodes={'4.1.2.60'};
rxnsToAdd.rxnComps=[1];
rxnsToAdd.lb=[-50];
rxnsToAdd.ub=[50];
rxnsToAdd.c=[0];
 
model13=addRxns(model12,rxnsToAdd,1);
[BBB, BBBprod]=checkBBBTasks(model13);
printRxnFormula(model13);
 
%   The met added previously (C21615) needs to be secreted
met_to_secrete={'C21615'}; %   Glycolaldehyde triphosphate
[model14 addedRxns]=addExchangeRxns(model13, 'OUT', met_to_secrete);
[BBB, BBBprod]=checkBBBTasks(model14);
 
%   Addition of the reaction R03503 (2.7.6.3)
rxnsToAdd.rxns={'R03503'};
rxnsToAdd.equations={'C00002 + C01300 <=> C00020 + C04807'}
rxnsToAdd.rxnNames={'ATP:6-hydroxymethyl-7,8-dihydropterin 6'};
rxnsToAdd.eccodes={'2.7.6.3'};
rxnsToAdd.rxnComps=[1];

rxnsToAdd.lb=[-50];
rxnsToAdd.ub=[50];
rxnsToAdd.c=[0];
 
model15=addRxns(model14,rxnsToAdd,1);
[BBB, BBBprod]=checkBBBTasks(model15);
 
%   In the One carbon pool by folate, adding of the following reaction
rxnsToAdd.rxns={'R00939'};
rxnsToAdd.equations={'C00101 + C00006 <=> C00415 + C00005 + C00080'}
rxnsToAdd.rxnNames={'5,6,7,8-tetrahydrofolate:NADP+ oxidoreductase'};
rxnsToAdd.eccodes={'1.5.1.3'};
rxnsToAdd.rxnComps=[1];
rxnsToAdd.lb=[-50];
rxnsToAdd.ub=[50];
rxnsToAdd.c=[0];

model16=addRxns(model15,rxnsToAdd,1);
[BBB, BBBprod]=checkBBBTasks(model16);

%   Addition of the reaction to produce the 10-Formyltetrahydrofolate
rxnsToAdd.rxns={'R00943'};
rxnsToAdd.equations={'C00101 + C00058 + C00002 <=> C00008 + C00009 + C00234'}
rxnsToAdd.rxnNames={'Formate:tetrahydrofolate ligase (ADP-forming)'};
rxnsToAdd.eccodes={'6.3.4.3'};
rxnsToAdd.rxnComps=[1];
rxnsToAdd.lb=[-50];
rxnsToAdd.ub=[50];
rxnsToAdd.c=[0];
 
model17=addRxns(model16,rxnsToAdd,1);
[BBB, BBBprod]=checkBBBTasks(model17);

%   Uptake and Secretion to produce methionine and S-Adenosyl-L-methionine
mets_to_uptake = {'C04489'}; %   5-Methyltetrahydropteroyltri-L-glutamate
[model18 addedRxns]=addExchangeRxns(model17,'BOTH',mets_to_uptake);

mets_to_secrete = {'C04144'}; % Tetrahydropteroyltri-L-glutamate
[model19 addedRxns]=addExchangeRxns(model18,'BOTH',mets_to_secrete);

%   Uptake C00229 (acyl-carrier protein)
metsToUptake = {'C00229'};
[model20 addedRxns]=addExchangeRxns(model19, 'BOTH',metsToUptake); 
[BBB, BBBprod]=checkBBBTasks(model20);

%   Secretion of Biomass_c to output secretion from the model
[model21 addedRxns]=addExchangeRxns(model20,'OUT','Biomass_c');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Growth in function of glucose and oxygen uptake rate
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   Varying the glucose uptake to test the effect of its variation on the growth of the cell
model21.lb(find(ismember(model21.rxns,'EXC_BOTH_C00031')))=0
sol_FBA=optimizeCbModel(model21); %   Solve FBA of the new model to see if the bacteria is able to grow

model21.lb(find(ismember(model21.rxns,'EXC_BOTH_C00031')))=-1
sol_FBA=optimizeCbModel(model21); %   Solve FBA of the new model to see if the bacteria is able to grow

model21.lb(find(ismember(model21.rxns,'EXC_BOTH_C00031')))=-5
sol_FBA=optimizeCbModel(model21); %   Solve FBA of the new model to see if the bacteria is able to grow

model21.lb(find(ismember(model21.rxns,'EXC_BOTH_C00031')))=-10
sol_FBA=optimizeCbModel(model21); %   Solve FBA of the new model to see if the bacteria is able to grow

model21.lb(find(ismember(model21.rxns,'EXC_BOTH_C00031')))=-15
sol_FBA=optimizeCbModel(model21); %   Solve FBA of the new model to see if the bacteria is able to grow

model21.lb(find(ismember(model21.rxns,'EXC_BOTH_C00031')))=-20
sol_FBA=optimizeCbModel(model21); %   Solve FBA of the new model to see if the bacteria is able to grow

%   Varying the oxygen uptake to test the effect of its variation on the growth of the cell
model21.lb(find(ismember(model21.rxns,'EXC_BOTH_C00007')))=-50
sol_FBA=optimizeCbModel(model21); %   Solve FBA of the new model to see if the bacteria is able to grow

model21.lb(find(ismember(model21.rxns,'EXC_BOTH_C00007')))=0
sol_FBA=optimizeCbModel(model21); %   Solve FBA of the new model to see if the bacteria is able to grow


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Flux variability analysis (FVA)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%   After having added the functions runMinMax.m and checkDirectionality.m 
%   to Raven, we ran the following code to perform the FVA
minmax = runMinMax(model21);
description_mm=checkDirectionality(model21,minmax(:,1),minmax(:,2),1,0);
