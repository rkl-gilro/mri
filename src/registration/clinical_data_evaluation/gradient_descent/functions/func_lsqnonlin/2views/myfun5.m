function F = myfun5(mesh0_X,in,out)

mesh1 = mesh0_X(1:length(mesh0_X)/2);
mesh2 = mesh0_X(length(mesh0_X)/2 + 1 :end);

mesh11 = reshape(mesh1',3,length(mesh0_X)/6);
mesh21 = reshape(mesh2',3,length(mesh0_X)/6);

tri1 = TriRep(in.source_tri.Triangulation,mesh11'); % define the new mesh for the axial mesh0_X(1:in.size_source,:)
tri2 = TriRep(in.source_tri.Triangulation,mesh21'); % define the new mesh for the sagittal mesh0_X(in.size_source+1:end,:)

% tri1 = TriRep(in.source_tri.Triangulation,mesh0_X(1:size(in.source_tri.X,1),:)); % define the new mesh for the axial
% tri2 = TriRep(in.source_tri.Triangulation,mesh0_X(size(in.source_tri.X,1)+1:end,:)); % define the new mesh for the sagittal


b2c_ncoord1 = baryToCart(tri1, in.current_tr, in.c2b_coord);
b2c_ncoord2 = baryToCart(tri2, in.current_tr, in.c2b_coord);

incr = in.n + in.size_source;

F = out.F;
sub_3 = in.sub_3;
sub_2 = in.sub_2;
axial_m1_par = in.axial_m1_par;
sag_m1_par = in.sag_m1_par;
vol_ax = in.vol_ax;
vol_sag = in.vol_sag;
rows = in.rows;
cols = in.cols;

parfor i = 1:in.n
    
    %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1 = axial_m1_par{i} * [b2c_ncoord1(i,:) 1]'; % 3D point to 2D point in the frame coordinates
    
    fl  = floor(tmp_v1(1) + 1);
    fl2 = floor(tmp_v1(2) + 1);
    cl  = ceil(tmp_v1(1) + 1);
    cl2 = ceil(tmp_v1(2) + 1);
    
    min_max_r = min(max(fl2,1),rows);
    min_max_r2 = min(max(cl2,1),rows);
    min_max_c = min(max(fl,1),cols);
    min_max_c2 = min(max(cl,1),cols);
    
    neig = [vol_ax(min_max_r,min_max_c,sub_3(i))   vol_ax(min_max_r,min_max_c2,sub_3(i));...
            vol_ax(min_max_r2,min_max_c,sub_3(i))  vol_ax(min_max_r2,min_max_c2,sub_3(i))];
    
    new_im_ax = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
    
    %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1 = sag_m1_par{i} * [b2c_ncoord2(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}
    
    fl  = floor(tmp_v1(1) + 1);
    fl2 = floor(tmp_v1(2) + 1);
    cl  = ceil(tmp_v1(1) + 1);
    cl2 = ceil(tmp_v1(2) + 1);

    min_max_r = min(max(fl2,1),rows);
    min_max_r2 = min(max(cl2,1),rows);
    min_max_c = min(max(fl,1),cols);
    min_max_c2 = min(max(cl,1),cols);
    
    neig = [vol_sag(min_max_r,min_max_c,sub_2(i))   vol_sag(min_max_r,min_max_c2,sub_2(i));...
            vol_sag(min_max_r2,min_max_c,sub_2(i))  vol_sag(min_max_r2, min_max_c2,sub_2(i))];
    
    new_im_sag = bilinear_interpolation(tmp_v1(2),tmp_v1(1),double(neig));
    
    
    %% Function
    F(i) = new_im_ax - new_im_sag;
    
end


% for i = 1:size(in.source_tri.X,1)
%     
%     F(i + in.n) = in.sqrt_lambda*(sum( abs(tri1.X(i,:) - mean(tri1.X(in.list_edges{i},:))) ));
%     F(i + incr) = in.sqrt_lambda*(sum( abs(tri2.X(i,:) - mean(tri2.X(in.list_edges{i},:))) ));
% end

% F = out.F;

