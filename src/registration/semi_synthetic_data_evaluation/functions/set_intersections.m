global optimizer

global t
t = 0:1/optimizer.t:1;

global vol_ax_eval
global vol_sag_eval
global vol_cor_eval

global X_ax_def
global Y_ax_def
global Z_ax_def

global X_sag_def
global Y_sag_def
global Z_sag_def

global X_cor_def
global Y_cor_def
global Z_cor_def

disp('--------- Calculate the s1 x s2 intersections between planes of the 2 given directions -----')
%% Intersections Axial & Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global var_cell1_v
global var_array1_v

[var_cell1_v, var_array1_v] = calculate_intersections(X_ax_def, Y_ax_def, Z_ax_def, X_sag_def, Y_sag_def, Z_sag_def, t, size(vol_ax_eval,3), size(vol_sag_eval,3));

disp('--------- Calculate the s1 x s3 intersections between planes of the 2 given directions -----')
%% Intersections Axial & Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global var_cell2_v
global var_array2_v

[var_cell2_v, var_array2_v] = calculate_intersections(X_ax_def, Y_ax_def, Z_ax_def, X_cor_def, Y_cor_def, Z_cor_def, t, size(vol_ax_eval,3), size(vol_cor_eval,3));

disp('--------- Calculate the s2 x s3 intersections between planes of the 2 given directions -----')
%% Intersections Axial & Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


global var_cell3_v
global var_array3_v

options = optimset('Display','off');

[var_cell3_v, var_array3_v] = calculate_intersections(X_cor_def, Y_cor_def, Z_cor_def, X_sag_def, Y_sag_def, Z_sag_def, t, size(vol_cor_eval,3), size(vol_sag_eval,3));

disp('--------- Calculate the source control points ( # N^3 ) -----')
%% Calculate the bounding box %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% bb = [xmin xmax;ymin ymax;zmin zmax]

global var_array

var_array = [var_array1_v;var_array2_v;var_array3_v];

% bb = [min(min(var_array(:,1)))-120 max(max(var_array(:,1)))+120; ...
%       min(min(var_array(:,2)))-120 max(max(var_array(:,2)))+120; ...
%       min(min(var_array(:,3)))-120 max(max(var_array(:,3)))+120];

bb = [min([min(X_cor_def(:)), min(X_sag_def(:)), min(X_ax_def(:))])-15 max([max(X_cor_def(:)), max(X_sag_def(:)), max(X_ax_def(:))])+15;...
      min([min(Y_cor_def(:)), min(Y_sag_def(:)), min(Y_ax_def(:))])-15 max([max(Y_cor_def(:)), max(Y_sag_def(:)), max(Y_ax_def(:))])+15;...
      min([min(Z_cor_def(:)), min(Z_sag_def(:)), min(Z_ax_def(:))])-15 max([max(Z_cor_def(:)), max(Z_sag_def(:)), max(Z_ax_def(:))])+15];
  
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
% nx = nx + 2;
% ny = ny + 2;
% nz = nz + 2;
% 
% %% Inner 
% l_x_in = l_x;
% l_y_in = l_y;
% l_z_in = l_z;
% %% Outer
% l_x_out = [-abs(l_x(1) - l_x(2)) + l_x(1) l_x abs(l_x(1) - l_x(2)) + l_x(end)];
% l_y_out = [-abs(l_y(1) - l_y(2)) + l_y(1) l_y abs(l_y(1) - l_y(2)) + l_y(end)];
% l_z_out = [-abs(l_z(1) - l_z(2)) + l_z(1) l_z abs(l_z(1) - l_z(2)) + l_z(end)];
% 
% %% Define the new source of control points
% l_x = l_x_n;
% l_y = l_y_n;
% l_z = l_z_n;

global source_control
global out_ind

source_control = zeros(nx * ny * nz, 3);

out_ind = zeros(nx * ny * nz, 1);

for i = 1:nx
    for j = 1:ny
        
        tmp =  1:nz;
        
        s2ind =  tmp + nz*(j-1 + ny*(i-1));
        
        source_control(s2ind,1) = repmat(l_x(i),nz,1);
        source_control(s2ind,2) = repmat(l_y(j),nz,1);
        source_control(s2ind,3) = l_z(tmp);
        
    end
end

for i = 1:nx
    for j = 1:ny
        
        for tmp =  1:nz
        
            s2ind =  tmp + nz*(j-1 + ny*(i-1));
            
            
            if i==1 || i==nx || j==1 || j==ny || tmp==1 || tmp==nz
                
                out_ind(s2ind) = 1;
                
            end
        end
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
global tetra2
tetra2 = [];

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

for i=1:nx-3
    
    for j=1:ny-3
        
        for k=1:nz-3
            
            ind_tmp =  k + (nz-2)*(j-1 + (ny-2)*(i-1));
            
            ind_tmp2 =  k+1 + (nz-2)*(j-1 + (ny-2)*(i-1));
            ind_tmp3 =  k + (nz-2)*(j-1 + (ny-2)*(i));
            ind_tmp4 =  k + (nz-2)*(j + (ny-2)*(i));
            ind_tmp5 =  k + (nz-2)*(j + (ny-2)*(i-1));
            ind_tmp6 =  k+1 + (nz-2)*(j + (ny-2)*(i));
            ind_tmp7 =  k+1 + (nz-2)*(j + (ny-2)*(i-1));
            ind_tmp8 =  k+1 + (nz-2)*(j-1 + (ny-2)*(i));
            
            tetra2 = [tetra2;...
                    ind_tmp  ind_tmp2 ind_tmp3 ind_tmp5;...
                    ind_tmp2 ind_tmp3 ind_tmp5 ind_tmp7;...
                    ind_tmp3 ind_tmp7 ind_tmp8 ind_tmp2;...
                    ind_tmp4 ind_tmp7 ind_tmp8 ind_tmp3;...
                    ind_tmp5 ind_tmp7 ind_tmp3 ind_tmp4;...
                    ind_tmp6 ind_tmp7 ind_tmp8 ind_tmp4];

        end
        
    end
    
end

% alpha(.1)
% axis equal

disp('--------- Convert the computed mesh into triangulation class -----')
%% Convert the computed mesh into a triangulation class is gonna help us for computing the
%% vertices or tetrahedrons that contain a query of points

global source_tri_v
global source_tri_v2
global list_edges_v
global list_edges_v2

trep = TriRep(tetra, source_control);
source_tri_v = trep;

% trep = TriRep(tetra2, source_control(out_ind==0,:));
% source_tri_v2 = trep;

list_edges_v = edges_connected(source_tri_v);

% list_edges_v2 = edges_connected(TriRep(tetra2, source_control(out_ind==0,:)));

%% Calculate the barycentric coordinates for the intersection points 
disp('--------- Preparing barycentric coordinates for first intersection: axial/sagittal  -----')
size1_int1 = size(var_cell1_v,3);
size2_int1 = size(var_cell1_v,2);
size3_int1 = size(var_cell1_v,1);

global sub_1_int1_v
global sub_2_int1_v
global sub_3_int1_v
global current_tr_int1_v
global c2b_coord_int1_v

[sub_1_int1_v, sub_2_int1_v, sub_3_int1_v] = ind2sub([size1_int1 size2_int1 size3_int1],1:size(var_array1_v,1));
[current_tr_int1_v, c2b_coord_int1_v]    = tsearchn(source_tri_v.X, source_tri_v.Triangulation, [var_array1_v(:,1) var_array1_v(:,2) var_array1_v(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array

disp('--------- Preparing barycentric coordinates for second intersection: axial/coronal  -----')
size1_int2 = size(var_cell2_v,3);
size2_int2 = size(var_cell2_v,2);
size3_int2 = size(var_cell2_v,1);

global sub_1_int2_v
global sub_2_int2_v
global sub_3_int2_v
global current_tr_int2_v
global c2b_coord_int2_v

[sub_1_int2_v, sub_2_int2_v, sub_3_int2_v] = ind2sub([size1_int2 size2_int2 size3_int2],1:size(var_array2_v,1));
[current_tr_int2_v, c2b_coord_int2_v]    = tsearchn(source_tri_v.X,source_tri_v.Triangulation,[var_array2_v(:,1) var_array2_v(:,2) var_array2_v(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
  

disp('--------- Preparing barycentric coordinates for third intersection: sagittal/coronal  -----')
size1_int3 = size(var_cell3_v,3);
size2_int3 = size(var_cell3_v,2);
size3_int3 = size(var_cell3_v,1);

global sub_1_int3_v
global sub_2_int3_v
global sub_3_int3_v
global current_tr_int3_v
global c2b_coord_int3_v


[sub_1_int3_v, sub_2_int3_v, sub_3_int3_v] = ind2sub([size1_int3 size2_int3 size3_int3],1:size(var_array3_v,1));
[current_tr_int3_v, c2b_coord_int3_v]    = tsearchn(source_tri_v.X,source_tri_v.Triangulation,[var_array3_v(:,1) var_array3_v(:,2) var_array3_v(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array

%% Calculate the barycentric coordinates for the intersection points 
% disp('--------- Preparing barycentric coordinates for first intersection: axial/sagittal  -----')
% size1_int1 = size(var_cell1_v,3);
% size2_int1 = size(var_cell1_v,2);
% size3_int1 = size(var_cell1_v,1);
% 
% global sub_1_int1_v2
% global sub_2_int1_v2
% global sub_3_int1_v2
% global current_tr_int1_v2
% global c2b_coord_int1_v2
% 
% [sub_1_int1_v2, sub_2_int1_v2, sub_3_int1_v2] = ind2sub([size1_int1 size2_int1 size3_int1],1:size(var_array1_v,1));
% [current_tr_int1_v2, c2b_coord_int1_v2]    = tsearchn(source_tri_v2.X, source_tri_v2.Triangulation, [var_array1_v(:,1) var_array1_v(:,2) var_array1_v(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
% 
% disp('--------- Preparing barycentric coordinates for second intersection: axial/coronal  -----')
% size1_int2 = size(var_cell2_v,3);
% size2_int2 = size(var_cell2_v,2);
% size3_int2 = size(var_cell2_v,1);
% 
% global sub_1_int2_v2
% global sub_2_int2_v2
% global sub_3_int2_v2
% global current_tr_int2_v2
% global c2b_coord_int2_v2
% 
% [sub_1_int2_v2, sub_2_int2_v2, sub_3_int2_v2] = ind2sub([size1_int2 size2_int2 size3_int2],1:size(var_array2_v,1));
% [current_tr_int2_v2, c2b_coord_int2_v2]    = tsearchn(source_tri_v2.X,source_tri_v2.Triangulation,[var_array2_v(:,1) var_array2_v(:,2) var_array2_v(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array
%   
% 
% disp('--------- Preparing barycentric coordinates for third intersection: sagittal/coronal  -----')
% size1_int3 = size(var_cell3_v,3);
% size2_int3 = size(var_cell3_v,2);
% size3_int3 = size(var_cell3_v,1);
% 
% global sub_1_int3_v2
% global sub_2_int3_v2
% global sub_3_int3_v2
% global current_tr_int3_v2
% global c2b_coord_int3_v2
% 
% 
% [sub_1_int3_v2, sub_2_int3_v2, sub_3_int3_v2] = ind2sub([size1_int3 size2_int3 size3_int3],1:size(var_array3_v,1));
% [current_tr_int3_v2, c2b_coord_int3_v2]    = tsearchn(source_tri_v2.X,source_tri_v2.Triangulation,[var_array3_v(:,1) var_array3_v(:,2) var_array3_v(:,3)]); % in our case it is a scalar, to make it general we consider current_tr as an array


