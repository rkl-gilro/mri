%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Main program for T2 axial, sagittal and coronal alignment
%%
%% NOTE: change line 13 & 18 to change the folder and file names
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% clc;
close all;

global save_folder

save_folder = 'resources/registration/dallot/56gradt10_74_.01';

mkdir(save_folder);

%% Change the name of the file for saving the results
save_name = strcat(save_folder, '/dallot.mat');

%% Set up the parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-----------------------------------')
disp('--------- Set up parameters -------')
disp('-----------------------------------')

global optimizer
global t_reg
global measure_errors

t_reg = 0;
measure_errors = [];

%% Set up the number of points per intersectin 't', 
%% 'lambda' for the regularization term
%% 'nx, ny, nz' for the grid dimensions
optimizer.t     = 74;
optimizer.lambda = .01;
optimizer.nxyz  = [7 7 7];

if ~exist('set_up','var')
    set_parameters;
end

%% Start the optimization %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-----------------------------------')
disp('--------- Start Optimization  -----')
disp('-----------------------------------')

%% Using gradient descent
optimizer.steepestdescent = 1;
optimizer.grad       = 1; %% 0/1 without gradient/provide gradient
optimizer.maxiter    = 40; %% maximum of iterations for the optimzer
optimizer.checkderiv = 'off'; %% 'on'/'off' check or not the provided gradient with MATLAB one
optimizer.tolfun     = 0.1; 
optimizer.tolx       = 0.1; 
optimizer.con        = 0; %% 0/1 -> fminunc / fmincon

%% Using metaheuristics, differential evolution
optimizer.metaheuristic = 0;

%% Set the number of views (default 3)
optimizer.views         = 3; %% provide the number of views used for the optimization, 2 or 3

start_optimization( );

%% Calculate the new images from the optimization result

disp('--------- Compute the new images and visualization --------------')
compute_new_images( );

%% Calculate the measurements
disp('--------- Compute error measurements --------------')
measurements(  );

%% Save the output if the name is provided
save(save_name, 'source_tri', 'target_tri_ax' , 'target_tri_sag', 'target_tri_cor',...
                'axial_m','axial_m1','sag_m','sag_m1','cor_m','cor_m1', ...
                'new_axial', 'new_sagittal', 'new_coronal', ...
                'vol_ax', 'vol_sag' , 'vol_cor',...
                'var_array1', 'var_array2', 'var_array3',...
                'var_cell1', 'var_cell2', 'var_cell3', ...
                'optimizer', 't_reg', 'measure_errors');
            
%% Choose the images to visualize, i.e, [9 8 8] corresponds to slice 9 of the axial view,
%% slice 8 from sagittal view and slice 8 from coronal one.
show_slices = [19 9 12];
visualization( views, show_slices );

%% Clear workspace            
clearvars -except views dcmdir
