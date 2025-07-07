%% Installation of all necessary toolboxes (optPBN/BNPBN/SBToolbox2[optimisation))

% Assign current directory
fileposition = which(mfilename);
optpbn_directory = regexprep(fileposition, [mfilename '.m'],'');
current_directory=optpbn_directory;

%% 1) Add optimisation toolbox from SBToolbox2 to Matlab path

% Assign scripts directory and add to Matlab path
SBTB2_directory=[optpbn_directory filesep 'SBTOOLBOX2' filesep];
addpath(genpath(SBTB2_directory));
disp('Optimisation toolbox of SBToolbox2 is added to Matlab path')

%% 2) Add BNPBN toolbox to Matlab path

% Assign scripts directory and add to Matlab path
BNPBN_directory=[optpbn_directory filesep 'pbn-matlab-toolbox' filesep];
addpath(BNPBN_directory);
disp('BNPBN toolbox is added to Matlab path')

%% 3) Add optPBN toolbox to Matlab path

% Assign scripts directory and add to Matlab path
optPBN_directory=[optpbn_directory filesep 'optPBN_toolbox' filesep];
addpath(optPBN_directory)
disp('optPBN scripts are added to Matlab path')

%% Clear all intermediate variables and savepath
clear all
savepath
disp('The installation is done!')
