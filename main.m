clear; 
clc;
curr_path = pwd;
addpath(genpath(curr_path));
%% Experiment parameter %%
parameters;
%% TRILC (without trial reduction)
run_TRILC;
%% TRILC (with trial reduction)
run_TRILC_reduction;
%% data-driven optimal ILC (DDOILC)
run_DDOILC;
%% Two points simultaneous perturbation stochastic approximation (SPSA2)
run_SPSA2;
%% One point simultaneous perturbation stochastic approximation (SPSA1)
run_SPSA1;
%% Global minimum using inverse distance weighting and surrogate radial basis functions (GLIS) 
% % % % % % % % % % % % % % % % % % % % % % % % 
% GLIS requires to install additional package, please refer to the
% cited paper and follow their instructions. 
% Use our function: calfun.m, and then run GLIS main function.
% % % % % % % % % % % % % % % % % % % % % % % %
%% BOBYQA 
% % % % % % % % % % % % % % % % % % % % % % % % 
% BOBYQA requires to install additional package, please refer to the 
% cited papers and follow their instructions. 
% Use our function: calfun.m, and then run BOBYQA main function.
% % % % % % % % % % % % % % % % % % % % % % % % 
%% Visualization
picture;
