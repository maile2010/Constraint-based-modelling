%Load data
model = readCbModel('D:/UM/MATLAB/e_coli_core.mat')

metabolite_1 = model.mets{26};  
metabolite_2 = model.mets{55};  
reaction = model.rxns{4};
printRxnFormula(model, 'rxnAbbrList', reaction, 'printFlag', true); % name of reaction cores column 4
fprintf('Metabolite name: %s\n', metabolite_1); % name of meta 1 cores row 26
fprintf('Metabolite name: %s\n', metabolite_2); % name of meta 2 cores row 55
stoich_coeff = model.S(26, 4); 
%gene associated to reaction 4 row = reaction, colume = gen, 1 = gene is
%involved in this reaction, 0 is not
genes_associated = find(e_coli_core.rxnGeneMat(4, :) == 1);
genes = e_coli_core.genes(genes_associated);
disp(genes);

%% FBA
% Identify the reaction ID for glucose exchange
rxnID = 'EX_glc__D_e'; 

% Set the lower bound for glucose uptake
model.lb(strcmp(model.rxns, rxnID)) = -20; % Set to 20 mmol/gDW/h

% Set the objective function to biomass reaction
model.c = zeros(length(model.rxns), 1);  % Set all objective coefficients to 0
model.c(strcmp(model.rxns, 'BIOMASS_Ecoli_core_w_GAM')) = 1; % Set biomass as objective

% Solve the FBA problem
model_opt = optimizeCbModel(model);

% Calculate optimal growth rate
optimalGrowthRate = model_opt.f; 

%Anaerobic condition
%Set lower bound and  upper bound = 0
rxnID_2 = 'EX_o2_e'; 
model.lb(strcmp(model.rxns, rxnID_2)) = 0; % Set to 0
model.ub(strcmp(model.rxns, rxnID_2)) = 0; % Set to 0
% Set the objective function to biomass reaction
model.c = zeros(length(model.rxns), 1); 
model.c(strcmp(model.rxns, 'BIOMASS_Ecoli_core_w_GAM')) = 1;
% Solve the FBA problem
model_opt_2 = optimizeCbModel(model);
% Calculate optimal growth rate
optimalGrowthRate_2 = model_opt_2.f; 

%% half-aerobic condition
%Set lower bound and  upper bound = 0
% Set the upper bound for glucose uptake
rxnID = 'EX_glc__D_e'; 
model.lb(strcmp(model.rxns, rxnID)) = -20; % Set to -20 mmol/gDW/h
rxnID_2 = 'EX_o2_e'; 
model.lb(strcmp(model.rxns, rxnID_2)) = -41.75/2;
model.ub(strcmp(model.rxns, rxnID_2)) = 1000; 
% Set the objective function to biomass reaction
model.c = zeros(length(model.rxns), 1); 
model.c(strcmp(model.rxns, 'BIOMASS_Ecoli_core_w_GAM')) = 1;
% Solve the FBA problem
model_opt_3 = optimizeCbModel(model);
% Calculate optimal growth rate
optimalGrowthRate_3 = model_opt_3.f; 

%% maximize ATP 
%aerobic condition
% Set the upper bound for glucose uptake
rxnID = 'EX_glc__D_e'; 
model.lb(strcmp(model.rxns, rxnID)) = -20; % Set to 20 mmol/gDW/h
rxnID_2 = 'EX_o2_e'; 
model.lb(strcmp(model.rxns, rxnID_2)) = -1000; % Set to 0
model.ub(strcmp(model.rxns, rxnID_2)) = 1000; % Set to 0
% Set the objective function to biomass reaction
model.c = zeros(length(model.rxns), 1); 
model.c(strcmp(model.rxns, 'ATPM')) = 1;
% Solve the FBA problem
model_opt_4 = optimizeCbModel(model);
% Calculate optimal growth rate
optimalATP_1 = model_opt_4.f; 

%annerobic
rxnID = 'EX_glc__D_e'; 
model.lb(strcmp(model.rxns, rxnID)) = -20; % Set to 20 mmol/gDW/h
rxnID_2 = 'EX_o2_e'; 
model.lb(strcmp(model.rxns, rxnID_2)) = -0; % Set to 0
model.ub(strcmp(model.rxns, rxnID_2)) = 0; % Set to 0
% Set the objective function to biomass reaction
model.c = zeros(length(model.rxns), 1); 
model.c(strcmp(model.rxns, 'ATPM')) = 1;
% Solve the FBA problem
model_opt_5 = optimizeCbModel(model);
% Calculate optimal growth rate
optimalATP_2 = model_opt_5.f; 
