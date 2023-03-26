
%Generate KAB-XylABE strain model
%Writes new modified .mat model (iJN1463_KAB.mat) from P. putida GEM model iJN1463

%Load and prepare GEM model
modelBigg = 'iJN1463.mat'; 
modelbase=readCbModel(modelBigg);
modelkab = modelbase;

%Confirm Gurobi model
solverName = 'gurobi';
solverType = 'LP';
changeCobraSolver(solverName, solverType);


%Removing empty reactions/metabolites
modelkab = removeTrivialStoichiometry(modelkab);


%Add TPA transport and conversion to PCA from tph operon:
modelkab = addReaction(modelkab, 'EX_tpa_e', 'metaboliteList', {'tpa_e'} ,...
'stoichCoeffList', -1, 'reversible', true); %TPA exchange

modelkab = addReaction(modelkab, 'TPA2p','metaboliteList', {'tpa_e', 'tpa_p' },...
'stoichCoeffList', [-1; 1], 'reversible', true); %TPA to periplasm

modelkab = addReaction(modelkab, 'TPA2c','metaboliteList', {'h_p', 'tpa_p', 'h_c', 'tpa_c' },...
'stoichCoeffList', [-1; -1; 1; 1], 'reversible', false); %TPA periplasm to cytosol

modelkab = addReaction(modelkab, 'TPHA','metaboliteList', {'tpa_c', 'nadph_c', 'o2_c','h_c', 'dcd_c', 'nadp_c'},...
'stoichCoeffList', [-1; -1; -1; -1; 1; 1], 'reversible', false); %TPA to DCD by tphA

modelkab = addReaction(modelkab, 'TPHB','metaboliteList', {'dcd_c', 'nadp_c', '34dhbz_c', 'co2_c', 'nadph_c', 'h_c'},...
'stoichCoeffList', [-1; -1; 1; 1; 1; 1], 'reversible', false); %DCD to PCA by tphB


%Original model contains open reactions for EG metabolism. Here we remove
%this pathway as to simulate its repression:
modelkab = removeRxns(modelkab, { 'GLXCL', 'HPYRI', 'TRSARr', 'HPYRRx', 'HPYRRy', 'GLYCK', 'GLYCLTDy','GLYCLTDx'});


%Add Xylose transport and assimilation reactions:
modelkab = addReaction(modelkab, 'EX_xyl_e', 'metaboliteList', {'xyl_e'} ,...
'stoichCoeffList', -1, 'reversible', true); 

modelkab = addReaction(modelkab, 'Xyl2p','metaboliteList', {'xyl_e', 'xyl_p' },...
'stoichCoeffList', [-1; 1], 'reversible', true); 

modelkab = addReaction(modelkab, 'Xyl2c','metaboliteList', {'h_p', 'xyl_p', 'h_c', 'xyl_c' },...
'stoichCoeffList', [-1; -1; 1; 1], 'reversible', false); 

modelkab = addReaction(modelkab, 'XYLA','metaboliteList', {'xyl_c', 'xylu_c'},...
'stoichCoeffList', [-1; 1], 'reversible', false);

modelkab = addReaction(modelkab, 'XYLB','metaboliteList', {'xylu_c', 'atp_c', 'xu5p__D_c', 'adp_c', 'h_c'},...
'stoichCoeffList', [-1; -1; 1; 1; 1], 'reversible', false); 

%Remove knocked-out glucose dehydrogenase reaction (gcd), redundant gluconate reactions:
modelkab = removeRxns(modelkab, {'GLCDpp', 'GAD2ktpp', 'BGLApp'});


%Reaction L-tryptophan decarboxylase (LTDCLNo1) was removed as it was shown to be non-functional (Koyanagi et al., 2012)
modelkab = removeRxns(modelkab, {'LTDCLNo1'}); 

% Product reactions:
% PHA formation reaction, combined PHAs:
modelkab = addReaction(modelkab, 'PHA',...
'metaboliteList', {'C80aPHA_c', 'C100aPHA_c', 'C120aPHA_c', 'C121aPHA_c', 'C140aPHA_c', 'mclPHA_c'},...
'stoichCoeffList', [-1; -1; -1; -1; -1; 1], 'reversible', true);

modelkab = addSinkReactions(modelkab, 'mclPHA_c');

%Recombinant insulin production reaction
modelkab = addReaction(modelkab, 'RecINS', ...
'metaboliteList', {'ala__L_c', 'arg__L_c', 'asn__L_c', 'asp__L_c', 'cys__L_c', 'gln__L_c', ...
'glu__L_c', 'gly_c', 'his__L_c', 'ile__L_c', 'leu__L_c', 'lys__L_c', 'met__L_c' ,'phe__L_c' , ...
'pro__L_c' ,'ser__L_c' ,'thr__L_c' ,'trp__L_c' ,'tyr__L_c' ,'val__L_c', 'atp_c', 'gtp_c','INS_c', 'adp_c', 'gdp_c', 'pi_c'},...
'stoichCoeffList', [-18; -4; -7; -8; -8; -6; -13; -20; -11; -7; -19; -11; -3; -7; -5; -10; -9; -2; -7; -13; -404.77; -404.77; 1; 404.77; 404.77; 809.53], 'reversible', false); 

%Add sink reaction to allow INS to leave the system 
modelkab = addSinkReactions(modelkab, 'INS_c');


%Output file:

%Write new mat model for KAB-XylABE, can change to desired format
%‘text’,’xls’, or ‘sbml’:

writeCbModel(modelkab, 'fileName', 'iJN1463_KAB.mat','format','mat')