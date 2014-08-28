

disp('--------- Convert the computed mesh into triangulation class -----')
%% Convert the computed mesh into a triangulation class is gonna help us for computing the
%% vertices or tetrahedrons that contain a query of points
global source_control
global source_tri
global list_edges

        
trep = TriRep(tetra, source_control);
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

