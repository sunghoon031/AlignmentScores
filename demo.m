folder = fileparts(which(mfilename)); 
addpath(genpath(folder));

clear all; close all; clc;


%% Notations:
% 1. xyz is the camera position viewed from the external reference frame.
% 2. R is the 3 x 3 matrix and its each column corresponds to the unit 
%    vector in x, y, z direction, viewed from the external reference frame.
%    This means that if the camera was rotated by some rotation M, its R
%    would be updated as M*R.

n_total = 100; % Number of poses to be evaluated
n_outlier = 30; % Number of outlier poses
sigma_xyz = 0.03; % Noise level in the camera translation estimation
sigma_R = 3; % Noise level in the camera rotation estimation (in deg)

[xyz_gt, R_gt, xyz_input, R_input] = GenerateSyntheticData(n_total, n_outlier, sigma_xyz, sigma_R);


%% Compute TAS, RAS, and PAS:

TAS = ComputeTAS(xyz_gt,xyz_input);

RAS = ComputeRAS(R_input, R_gt);

PAS = 0.5*(TAS+RAS);
 
disp(['TAS = ', num2str(TAS), ' , RAS = ', num2str(RAS) , ', PAS = ', num2str(PAS)])



