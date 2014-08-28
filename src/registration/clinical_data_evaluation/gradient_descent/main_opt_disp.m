clc;
close all;

global save_folder
global optimizer
global t_reg
global measure_errors
global source_tri
global source_control

t_reg = 0;
measure_errors = [];

%% Set up the parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
disp('-----------------------------------')
disp('--------- Set up parameters -------')
disp('-----------------------------------')

%% Set up the number of points per intersectin 't',
%% 'lambda' for the regularization term
%% 'nx, ny, nz' for the grid dimensions
optimizer.t      = 54;
optimizer.lambda = .01;
optimizer.nxyz   = 3.*ones(1,3);

set_parameters2;

% Normal distribution for x,y,z
mesh_size = size(source_control,1);

n = 1000;
a = -5;
b =  5;
cov = .5; % covariance
m   = 0; % mean
Xrand = a + (b-a).*rand(mesh_size,n) + m.*ones(mesh_size,n);
Yrand = a + (b-a).*rand(mesh_size,n) + m.*ones(mesh_size,n);
Zrand = a + (b-a).*rand(mesh_size,n) + m.*ones(mesh_size,n);



for iter = 4:n
        
        disp('-------------------------------------------------------------------------')
        disp(['----------- Iteration: ',num2str(iter),' ----------------------------------'])
        disp('-------------------------------------------------------------------------')
        disp(' ')
        

        %% Create folder for the new registration parameters
        save_folder = strcat('resources/registration/lassalas/random_disp/',num2str(iter),'t',num2str(optimizer.t),'_lambda',num2str(optimizer.lambda),...
            '_nxyz',num2str(optimizer.nxyz(1)), num2str(optimizer.nxyz(2)), num2str(optimizer.nxyz(3)));
        mkdir(save_folder);
        
        %% Change the name of the file for saving the results
        save_name = strcat(save_folder,'/lassalas.mat');
        
        set_parameters2;
        
        source_control(:,1) = source_control(:,1) + Xrand(:,iter);
        source_control(:,2) = source_control(:,2) + Yrand(:,iter);
        source_control(:,3) = source_control(:,3) + Zrand(:,iter);
        
        x_dis = Xrand(:,iter);
        y_dis = Yrand(:,iter);
        z_dis = Zrand(:,iter);
        
        %% Calculate the coordinates of the intersection points given the mesh
        set_barcoordinates( );
        
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
        optimizer.tolx       = 0.01;
        optimizer.con        = 0; %% 0/1 -> fminunc / fmincon
        
        %% Using metaheuristics, differential evolution
        optimizer.metaheuristic = 0;
        
        %% Set the number of views (default 3)
        optimizer.views         = 3; %% provide the number of views used for the optimization, 2 or 3
        
        start_optimization( );
        
        %% Calculate the new images from the optimization result
        
        disp('--------- Compute the new images and visualization --------------')
        compute_new_images( );

        %% Compute the error and measurements
        disp('--------- Compute error measurements ----------------------------')
        measurements( );
        
        %% Save the output if the name is provided
        save(save_name, 'source_tri', 'target_tri_ax' , 'target_tri_sag', 'target_tri_cor',...
                        'axial_m','axial_m1','sag_m','sag_m1','cor_m','cor_m1', ...
                        'new_axial', 'new_sagittal', 'new_coronal', ...
                        'vol_ax', 'vol_sag' , 'vol_cor',...
                        'optimizer', 't_reg', 'measure_errors', ...
                        'x_dis', 'y_dis', 'z_dis');
        
        %% Clear workspace
        clearvars -except save_folder optimizer measure_errors t_reg vol_ax vol_sag vol_cor ...
                          views dcmdir iter source_control Xrand Yrand Zrand x_dis y_dis z_dis
        source_control = [];

end
