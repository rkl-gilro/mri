clc;
close all;

global save_folder
global optimizer
global t_reg
global measure_errors

global target_tri_ax
global target_tri_sag
global target_tri_cor
global source_tri

t_reg = 0;
measure_errors = [];

% Number of iterations
n = 100000;

% Values of the function in each iteration
F = zeros(1, n);

parametrization = 54;
lambda = .01;
nxyz = 5;

%% Set up the parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-----------------------------------')
disp('--------- Set up parameters -------')
disp('-----------------------------------')

%% Set up the number of points per intersectin 't',
%% 'lambda' for the regularization term
%% 'nx, ny, nz' for the grid dimensions
optimizer.t      = parametrization;
optimizer.lambda = lambda;
optimizer.nxyz   = nxyz.*ones(1,3);

%% Create folder for the new registration parameters
save_folder = strcat('resources/registration/lassalas/montecarlo/2t',num2str(optimizer.t),'_lambda',num2str(optimizer.lambda),...
    '_nxyz',num2str(optimizer.nxyz(1)), num2str(optimizer.nxyz(2)), num2str(optimizer.nxyz(3)));
mkdir(save_folder);

if ~exist('set_up','var')
    set_parameters;
    set_up = 1;
end

% Normal distribution for x,y,z
mesh_size = size(source_tri.X,1);
a = -5;
b =  5;
cov = .5; % covariance
m   = 0; % mean
Xrand = a + (b-a).*rand(mesh_size,n) + m.*ones(mesh_size,n);
Yrand = a + (b-a).*rand(mesh_size,n) + m.*ones(mesh_size,n);
Zrand = a + (b-a).*rand(mesh_size,n) + m.*ones(mesh_size,n);

%% Set up the variables for energy functional
[~, ~] = preparing4;

for iter = 408:n
        
        disp('-------------------------------------------------------------------------')
        disp(['----------- Iteration: ',num2str(iter),' ----------------------------------'])
        disp('-------------------------------------------------------------------------')
        disp(' ')
        %% Compute the functional value %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        F(iter) = montecarlo(Xrand(:,iter), Yrand(:,iter), Zrand(:,iter));
        
        disp('-------------------------------------------------------------------------')
        disp(['Functional Value: ',num2str(F(iter))])
        disp('-------------------------------------------------------------------------')
        %% Calculate the new images from the optimization result
        
        disp('--------- Compute the new images and visualization --------------')
%         compute_new_images( );
        
        %% Choose the images to visualize, i.e, [9 8 8] corresponds to slice 9 of the axial view,
        %% slice 8 from sagittal view and slice 8 from coronal one.
        % show_slices = [19 9 12];
        % visualization( views, show_slices );
        
        %% Compute the error and measurements
        disp('--------- Compute error measurements ----------------------------')
        measurements( );
        
        %% Change the name of the file for saving the results
        [~, zeros_str] = add_zeros(iter, floor(log10( n )) + 1);
        save_name = strcat(save_folder, '/lassalas', zeros_str, num2str(iter), '.mat');
        
        X = Xrand(:,iter);
        Y = Yrand(:,iter);
        Z = Zrand(:,iter);
        
        %% Save the output if the name is provided
        save(save_name, 'source_tri', 'target_tri_ax' , 'target_tri_sag', 'target_tri_cor',...
                        'axial_m','axial_m1','sag_m','sag_m1','cor_m','cor_m1', ... %'new_axial', 'new_sagittal', 'new_coronal',..
                        'vol_ax', 'vol_sag' , 'vol_cor',...
                        'optimizer', 'iter', 'measure_errors', 'F', 'X', 'Y', 'Z');
        
        %% Clear workspace
        clearvars -except save_folder set_up optimizer measure_errors t_reg vol_ax vol_sag vol_cor views dcmdir p p_length  iter ...
                          Xrand Yrand Zrand n F cov m parametrization lambda nxyz source_tri target_tri_ax target_tri_sag target_tri_cor ...
                          axial_m axial_m1 sag_m sag_m1 cor_m cor_m1 ... %'new_axial', 'new_sagittal', 'new_coronal',..
                          vol_ax vol_sag vol_cor
    
end
