%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Imaging we have the deformation (target control points). 
%% Question: how to determine the location of a point p after the deformation,
%% p is not a source control point.
%% Answer: with barycentric coordinates
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%clc;
%close all;

global source_tri
global var_array1
global axial_m1
global sag_m1
global var_cell1
global vol_ax
global vol_sag

rows = size(vol_ax,1);
cols = size(vol_ax,2);

tri = source_tri.Triangulation;

%% Plot the original triangulation and reference points.
figure
subplot(1,2,1);
tetramesh(source_tri); hold on;
alpha(.1)
plot3(var_array1(:,1), var_array1(:,2),var_array1(:,3),'*r'); 
hold off;
axis equal;

%% Stretch the triangulation and compute the mapped locations of the incenters on the deformed triangulation.
% noise_data = source_tri.X + randn(size(source_tri.X));

% target_tri = TriRep(tri,source_tri.X(:,1), 3.*source_tri.X(:,2), 10 + source_tri.X(:,3));
% target_tri = TriRep(tri,noise_data(:,1),noise_data(:,2),noise_data(:,3));
target_tri  = TriRep(source_tri.Triangulation,xfinal(1:size(xfinal,1)/2,:));
target_tri2 = TriRep(source_tri.Triangulation,xfinal(size(xfinal,1)/2+1:end,:));

current_tr = tsearchn(source_tri.X,source_tri.Triangulation,[var_array1(:,1) var_array1(:,2) var_array1(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array

% New data positions with respect the new computed meshes
b = cartToBary(source_tri,current_tr,[var_array1(:,1) var_array1(:,2) var_array1(:,3)]); % barycentric coordinates of the points wrt source tri
xc = baryToCart(target_tri, current_tr, b); % cartesian coordinates of the points wrt target tri ( p' )
xc2 = baryToCart(target_tri2, current_tr, b);

%% Get the intensity values of the p' in the first direction
opt_im_ax = zeros(size(vol_ax));
opt_im_sag = zeros(size(vol_sag));

for i = 1:size(var_array1,1)
    
    [sub_tmp1 sub_tmp2 sub_tmp3] = ind2sub([size(var_cell1,3) size(var_cell1,2) size(var_cell1,1)],i);
    
    tmp_v1 = axial_m1{sub_tmp3} * [xc(i,:) 1]'; % 3D point to 2D point in the frame coordinates
    %tmp_v1 = axial_m1{sub_tmp3} * xc(i,:)';
    tmp_v = [tmp_v1(2) tmp_v1(1) tmp_v1(3)]; % switch the first and second coordinates, to have i',j' (rows,cols)
    
    %% Make sure the indexes are correct and use bilinear interpolation

    neig = [vol_ax(min(max(floor(tmp_v(1) + 1),1),rows),min(max(floor(tmp_v(2) + 1),1),cols),sub_tmp3) vol_ax(min(max(floor(tmp_v(1) + 1),1),rows),min(max(ceil(tmp_v(2) + 1),1),cols),sub_tmp3);...
            vol_ax(min(max(ceil(tmp_v(1) + 1),1),rows), min(max(floor(tmp_v(2) + 1),1),cols),sub_tmp3) vol_ax(min(max(ceil(tmp_v(1) + 1),1),rows), min(max(ceil(tmp_v(2) + 1),1),cols),sub_tmp3)];


    new_i = min(max(floor(tmp_v(1) + 1),1),rows);%round(tmp_v(1) + 1) % matlab starts always at 1
    new_j = min(max(floor(tmp_v(2) + 1),1),cols);%round(tmp_v(2) + 1) % matlab starts always at 1
    new_k = sub_tmp3;

    opt_im_ax(new_i,new_j,new_k) = bilinear_interpolation(tmp_v(1),tmp_v(2),double(neig));
    aa = opt_im_ax(new_i,new_j,new_k);
    %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tmp_v1 = sag_m1{sub_tmp2} * [xc2(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}
    %tmp_v1 = sag_m1{sub_tmp2} * xc2(i,:)';
    tmp_v = [tmp_v1(2) tmp_v1(1) tmp_v1(3)]; % switch the first and second coordinates, to have i',j' (rows,cols)
    
    %% Make sure the indexes are correct
    
    neig = [vol_sag(min(max(floor(tmp_v(1) + 1),1),rows),min(max(floor(tmp_v(2) + 1),1),cols),sub_tmp2) vol_sag(min(max(floor(tmp_v(1) + 1),1),rows),min(max(ceil(tmp_v(2) + 1),1),cols),sub_tmp2);...
            vol_sag(min(max(ceil(tmp_v(1) + 1),1),rows), min(max(floor(tmp_v(2) + 1),1),cols),sub_tmp2) vol_sag(min(max(ceil(tmp_v(1) + 1),1),rows), min(max(ceil(tmp_v(2) + 1),1),cols),sub_tmp2)];

    new_i = min(max(floor(tmp_v(1) + 1),1),rows); % matlab starts always at 1
    new_j = min(max(floor(tmp_v(2) + 1),1),cols); % matlab starts always at 1
    new_k = sub_tmp2;

    %% Make sure the indexes are correct

    opt_im_sag(new_i,new_j,new_k) = bilinear_interpolation(tmp_v(1),tmp_v(2),double(neig));
    s(i) = abs(opt_im_sag(new_i,new_j,new_k) - aa);
end 

mean(s)
std(s)
show_results(opt_im_ax);
show_results(opt_im_sag);
% Plot the deformed triangulation and mapped locations of the reference points.
subplot(1,2,2);
tetramesh(target_tri2);title('Second mesh');  
alpha(.1)
hold on;
plot3(xc2(:,1), xc2(:,2), xc2(:,3), '*r'); 
hold off;
axis equal;

%% Plot the original triangulation and reference points.
figure;
subplot(1,2,1);
tetramesh(source_tri); hold on;
alpha(.1)
plot3(var_array1(:,1), var_array1(:,2),var_array1(:,3),'*r'); 
hold off;
axis equal;

subplot(1,2,2);
tetramesh(target_tri);title('First mesh'); 
alpha(.1)
hold on;
plot3(xc(:,1), xc(:,2), xc(:,3), '*r'); 
hold off;
axis equal;

varlist = {'new_k','new_j','new_i','neig','sub_tmp1','sub_tmp2','sub_tmp3','tmp_v1','tmp_v','aa'};
clear(varlist{:})
clear varlist
