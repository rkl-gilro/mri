% clc

global set_up

set_up = 1;

global optimizer

global vol_ax
global vol_sag
global vol_cor

global new_axial
global new_sagittal
global new_coronal

global gradx_ax
global grady_ax
global gradx_sag
global grady_sag
global gradx_cor
global grady_cor

views.axial    = double(views.axial);
views.sagittal = double(views.sagittal);
views.coronal  = double(views.coronal);

new_axial    = zeros(size(views.axial));
new_sagittal = zeros(size(views.sagittal));
new_coronal  = zeros(size(views.coronal));

rows = size(views.axial,1);
cols = size(views.axial,2);

total_ax = size(views.axial,3);
total_sag = size(views.sagittal,3);
total_cor = size(views.coronal,3);

global plane_ax
global plane_sag
global plane_cor

global axial_m
global axial_m1

global t
t = 0:1/optimizer.t:1;

global ortho
ortho = 1;

disp('--------- Normalize values (future work: bias field correction) -----')
%% Normalize values
min_val = min([min(views.axial(:)) min(views.sagittal(:)) min(views.coronal(:))])
max_val = max([max(views.axial(:)) max(views.sagittal(:)) max(views.coronal(:))])

vol_ax  = views.axial;
vol_sag = views.sagittal;
vol_cor = views.coronal;

% vol_ax  = range(views.axial,    min_val, max_val);
% vol_sag = range(views.sagittal, min_val, max_val);
% vol_cor = range(views.coronal,  min_val, max_val);

disp('-----------------------------------')
disp('---------- Image denoising   ------')
disp('-----------------------------------')
%% Image denoising %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Rician noise
disp('-------------Rician Noise ---------')
M     = 3;
alpha = 1;
show  = 0;
h   = 16;

% if ~isempty('vol_ax')
%     [vol_ax, ~ ] = rician_noise(vol_ax,  [M, alpha, h], show );
%     [vol_sag, ~] = rician_noise(vol_sag, [M, alpha, h], show );
%     [vol_cor, ~] = rician_noise(vol_cor, [M, alpha, h], show );
% end

%% Anisotropic diffusion
% disp('------------------------------ Anisotropic Diffusion ---------')
% for i=1:total_ax
%     vol_ax(:,:,i) = anisodiff2D(views.axial(:,:,i), 10, 1/7, 30, 1);
% end
% 
% for i=1:total_sag
%     vol_sag(:,:,i) = anisodiff2D(views.sagittal(:,:,i), 10, 1/7, 30, 1);
% end
% for i=1:total_cor
%     vol_cor(:,:,i) = anisodiff2D(views.coronal(:,:,i), 10, 1/7, 30, 1);
% end

%% Volume bias correction
% disp('----------- Bias Correction --------')
% umbral = 1; % 0/1 don't erase background/erase background

% res_ax = [1 1 round(views.axial_info{1,1}.SpacingBetweenSlices / views.axial_info{1,1}.PixelSpacing(1)) ]; 
% [vol_ax, ~] = bias_correction_vol(vol_ax, res_ax, umbral);
% 
% res_sag = [1 1 round(views.sagittal_info{1,1}.SpacingBetweenSlices / views.sagittal_info{1,1}.PixelSpacing(1)) ]; 
% [vol_sag, ~] = bias_correction_vol(vol_sag, res_sag, umbral);
% 
% res_cor = [1 1 round(views.coronal_info{1,1}.SpacingBetweenSlices / views.coronal_info{1,1}.PixelSpacing(1)) ]; 
% [vol_cor, ~] = bias_correction_vol(vol_cor, res_cor, umbral);



% show_results(vol_ax)
% show_results(vol_sag)
% show_results(vol_cor)

disp('--------- Calculate image gradients -----')
%% Image gradients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = fspecial('prewitt');

% H = [  17    61    17  ;
%         0     0     0   ;
%         -17   -61   -17 ];
    
sigm = 5;
G = fspecial('gaussian',[2*sigm+1 2*sigm+1],sigm);
[dx, dy] = gradient(fspecial('gauss',[2*sigm+1 2*sigm+1],sigm)); % G is a 2D gaussain

        
for i=1:total_ax
    
    grady_ax(:,:,i) = imfilter(im2double(vol_ax(:,:,i)), dy,'same'); % conv2(vol_ax(:,:,i), H', 'same');%
    gradx_ax(:,:,i) = imfilter(im2double(vol_ax(:,:,i)), dx,'same');

end

for i=1:total_sag
    
    grady_sag(:,:,i) = imfilter(im2double(vol_sag(:,:,i)), dy, 'same');
    gradx_sag(:,:,i) = imfilter(im2double(vol_sag(:,:,i)), dx, 'same');

end
for i=1:total_cor
    
    grady_cor(:,:,i) = imfilter(im2double(vol_cor(:,:,i)), dy, 'same');
    gradx_cor(:,:,i) = imfilter(im2double(vol_cor(:,:,i)), dx, 'same');
    
end

disp('--------- Calculate M & M^(-1) and plane eq. for each slice in the first direction -----')
M = zeros(4,3);
%% Axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[axial_m axial_m1 X_ax Y_ax Z_ax plane_ax] = calculate_transforms(views.axial_info, size(views.axial),  [0 size(views.axial,2)-1], [0 size(views.axial,1)-1], ortho, 1);

disp('--------- Calculate  M & M^(-1) and plane eq. for each slice in the second direction -----')
%% Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global sag_m
global sag_m1

[sag_m sag_m1 X_sag Y_sag Z_sag plane_sag] = calculate_transforms(views.sagittal_info, size(views.sagittal), [0 size(views.sagittal,2)-1], [0 size(views.sagittal,1)-1], ortho, 2);

disp('--------- Calculate  M & M^(-1) and plane eq. for each slice in the third direction -----')
%% Coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global cor_m
global cor_m1

[cor_m cor_m1 X_cor Y_cor Z_cor plane_cor] = calculate_transforms(views.coronal_info, size(views.coronal), [0 size(views.coronal,2)-1], [0 size(views.coronal,1)-1],ortho, 3);

disp('--------- Calculate the s1 x s2 intersections between planes of the 2 given directions -----')
%% Intersections Axial & Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global var_cell1
global var_array1

[var_cell1, var_array1] = calculate_intersections(X_ax, Y_ax, Z_ax, X_sag, Y_sag, Z_sag, t, size(views.axial,3), size(views.sagittal,3));

disp('--------- Calculate the s1 x s3 intersections between planes of the 2 given directions -----')
%% Intersections Axial & Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global t2
t2 = 0:1/5:1;

global var_cell2
global var_array2

[var_cell2, var_array2] = calculate_intersections(X_ax, Y_ax, Z_ax, X_cor, Y_cor, Z_cor, t, size(views.axial,3), size(views.coronal,3));

disp('--------- Calculate the s2 x s3 intersections between planes of the 2 given directions -----')
%% Intersections Axial & Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global t3
t3 = 0:1/5:1;

global var_cell3
global var_array3

options = optimset('Display','off');

[var_cell3, var_array3] = calculate_intersections(X_cor, Y_cor, Z_cor, X_sag, Y_sag, Z_sag, t, size(views.coronal,3), size(views.sagittal,3));

disp('--------- Calculate the source control points ( # N^3 ) -----')
%% Calculate the bounding box %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bb = [xmin xmax;ymin ymax;zmin zmax]

global var_array

var_array = [var_array1;var_array2;var_array3];

% bb = [min(min(var_array(:,1)))-120 max(max(var_array(:,1)))+120; ...
%       min(min(var_array(:,2)))-120 max(max(var_array(:,2)))+120; ...
%       min(min(var_array(:,3)))-120 max(max(var_array(:,3)))+120];

bb = [min([min(X_cor(:)), min(X_sag(:)), min(X_ax(:))])-15 max([max(X_cor(:)), max(X_sag(:)), max(X_ax(:))])+15;...
      min([min(Y_cor(:)), min(Y_sag(:)), min(Y_ax(:))])-15 max([max(Y_cor(:)), max(Y_sag(:)), max(Y_ax(:))])+15;...
      min([min(Z_cor(:)), min(Z_sag(:)), min(Z_ax(:))])-15 max([max(Z_cor(:)), max(Z_sag(:)), max(Z_ax(:))])+15];
  
% Create the source control points
nx = optimizer.nxyz(1);
ny = optimizer.nxyz(2);
nz = optimizer.nxyz(3);

l_x = linspace(bb(1,1), bb(1,2), nx);
l_y = linspace(bb(2,1), bb(2,2), ny);
l_z = linspace(bb(3,1), bb(3,2), nz);

%% Add anothergrid (experimental)
% l_x_n = [-abs(l_x(1) - l_x(2)) + l_x(1) l_x abs(l_x(1) - l_x(2)) + l_x(end)];
% l_y_n = [-abs(l_y(1) - l_y(2)) + l_y(1) l_y abs(l_y(1) - l_y(2)) + l_y(end)];
% l_z_n = [-abs(l_z(1) - l_z(2)) + l_z(1) l_z abs(l_z(1) - l_z(2)) + l_z(end)];
% 
% nx = nx + 1;
% ny = nx + 1;
% nz = nx + 1;
% 
% l_x = l_x_n;
% l_y = l_y_n;
% l_z = l_z_n;

global source_control
source_control = zeros(nx * ny * nz, 3);

for i = 1:nx
    for j = 1:ny
        
        tmp =  1:nz;
        
        s2ind =  tmp + nz*(j-1 + ny*(i-1));
        
        source_control(s2ind,1) = repmat(l_x(i),nz,1);
        source_control(s2ind,2) = repmat(l_y(j),nz,1);
        source_control(s2ind,3) = l_z(tmp);
        
    end
end

% Plot the source control points
% for i = 1:nx * ny * nz
%     plot3(source_control(i,1),source_control(i,2),source_control(i,3),'k+');hold on
% end

disp('--------- Calculate the source control points mesh (tetrahedrons) -----')
%% Define the mesh for FEM %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global tetra
tetra = [];

% figure;

for i=1:nx-1
    
    for j=1:ny-1
        
        for k=1:nz-1
            
            ind_tmp =  k + nz*(j-1 + ny*(i-1));
            
            ind_tmp2 =  k+1 + nz*(j-1 + ny*(i-1));
            ind_tmp3 =  k + nz*(j-1 + ny*(i));
            ind_tmp4 =  k + nz*(j + ny*(i));
            ind_tmp5 =  k + nz*(j + ny*(i-1));
            ind_tmp6 =  k+1 + nz*(j + ny*(i));
            ind_tmp7 =  k+1 + nz*(j + ny*(i-1));
            ind_tmp8 =  k+1 + nz*(j-1 + ny*(i));
            
            tetra = [tetra;...
                    ind_tmp  ind_tmp2 ind_tmp3 ind_tmp5;...
                    ind_tmp2 ind_tmp3 ind_tmp5 ind_tmp7;...
                    ind_tmp3 ind_tmp7 ind_tmp8 ind_tmp2;...
                    ind_tmp4 ind_tmp7 ind_tmp8 ind_tmp3;...
                    ind_tmp5 ind_tmp7 ind_tmp3 ind_tmp4;...
                    ind_tmp6 ind_tmp7 ind_tmp8 ind_tmp4];
            
            % Just for the plotting
%             vari = [source_control(ind_tmp,:);source_control(ind_tmp2,:);source_control(ind_tmp3,:);source_control(ind_tmp4,:);...
%                     source_control(ind_tmp5,:);source_control(ind_tmp6,:);source_control(ind_tmp7,:);source_control(ind_tmp8,:)];
%             dt = DelaunayTri(vari);
%             
%             tetramesh(dt);
%             alpha(.1)
%             axis equal
%             axis off
        end
        
    end
    
end

% alpha(.1)
% axis equal

disp('--------- Convert the computed mesh into triangulation class -----')
%% Convert the computed mesh into a triangulation class is gonna help us for computing the
%% vertices or tetrahedrons that contain a query of points

global source_tri
global list_edges

trep = TriRep(tetra,source_control);
source_tri = trep;

list_edges = edges_connected(source_tri);

%% Calculate the barycentric coordinates for the intersection points 
disp('--------- Preparing barycentric coordinates for first intersection: axial/sagittal  -----')
size1_int1 = size(var_cell1,3);
size2_int1 = size(var_cell1,2);
size3_int1 = size(var_cell1,1);

global sub_1_int1
global sub_2_int1
global sub_3_int1
global current_tr_int1
global c2b_coord_int1

[sub_1_int1, sub_2_int1, sub_3_int1] = ind2sub([size1_int1 size2_int1 size3_int1],1:size(var_array1,1));
[current_tr_int1, c2b_coord_int1]    = tsearchn(source_tri.X, source_tri.Triangulation, [var_array1(:,1) var_array1(:,2) var_array1(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array

disp('--------- Preparing barycentric coordinates for second intersection: axial/coronal  -----')
size1_int2 = size(var_cell2,3);
size2_int2 = size(var_cell2,2);
size3_int2 = size(var_cell2,1);

global sub_1_int2
global sub_2_int2
global sub_3_int2
global current_tr_int2
global c2b_coord_int2

[sub_1_int2, sub_2_int2, sub_3_int2] = ind2sub([size1_int2 size2_int2 size3_int2],1:size(var_array2,1));
[current_tr_int2, c2b_coord_int2]    = tsearchn(source_tri.X,source_tri.Triangulation,[var_array2(:,1) var_array2(:,2) var_array2(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
  

disp('--------- Preparing barycentric coordinates for third intersection: sagittal/coronal  -----')
size1_int3 = size(var_cell3,3);
size2_int3 = size(var_cell3,2);
size3_int3 = size(var_cell3,1);

global sub_1_int3
global sub_2_int3
global sub_3_int3
global current_tr_int3
global c2b_coord_int3


[sub_1_int3, sub_2_int3, sub_3_int3] = ind2sub([size1_int3 size2_int3 size3_int3],1:size(var_array3,1));
[current_tr_int3, c2b_coord_int3]    = tsearchn(source_tri.X,source_tri.Triangulation,[var_array3(:,1) var_array3(:,2) var_array3(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array

%% Compute the arrays in case not all the intersections are going to be taken into account

global slice_ax
global slice_sag
global slice_cor

slice_ax  = 1:size(views.axial,3);
slice_sag = 1:size(views.sagittal,3);
slice_cor = 1:size(views.coronal,3);

global slice_int1
global slice_int2
global slice_int3

slice_int1 = [];
slice_int2 = [];
slice_int3 = [];

tmp_ind1 = 0;
tmp_ind2 = 1;
for i=1:length(slice_ax)
    for j=1:length(slice_sag)
        slice_int1(tmp_ind1+1:length(t)*tmp_ind2) = sub2ind([length(t), size(views.sagittal,3), size(views.axial,3)],   1:length(t), repmat(slice_sag(j),1,length(t)), repmat(slice_ax(i),1,length(t)));
        tmp_ind1 = length(t)*(tmp_ind2);
        tmp_ind2 = tmp_ind2 + 1;
    end
end


tmp_ind1 = 0;
tmp_ind2 = 1;
for i=1:length(slice_ax)
    for j=1:length(slice_cor)
        slice_int2(tmp_ind1+1:length(t)*tmp_ind2) = sub2ind([length(t), size(views.coronal,3),  size(views.axial,3)],   1:length(t), repmat(slice_cor(j),1,length(t)), repmat(slice_ax(i),1,length(t)));
        tmp_ind1 = length(t)*(tmp_ind2);
        tmp_ind2 = tmp_ind2 + 1;
    end
end

tmp_ind1 = 0;
tmp_ind2 = 1;
for i=1:length(slice_cor)
    for j=1:length(slice_sag)
        slice_int3(tmp_ind1+1:length(t)*tmp_ind2) = sub2ind([length(t), size(views.sagittal,3), size(views.coronal,3)], 1:length(t), repmat(slice_sag(j),1,length(t)), repmat(slice_cor(i),1,length(t)));
        tmp_ind1 = length(t)*(tmp_ind2);
        tmp_ind2 = tmp_ind2 + 1;
    end
end

