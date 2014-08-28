function [F, J] = myfun_ortho(mesh0_X)


global source_tri
global var_array1
global var_array2
global var_array3

global vol_ax
global vol_sag
global vol_cor

global sub_2_int1
global sub_3_int1

global sub_2_int2
global sub_3_int2

global sub_2_int3
global sub_3_int3
global current_tr_int1
global c2b_coord_int1
global current_tr_int2
global c2b_coord_int2
global current_tr_int3
global c2b_coord_int3
global list_edges
global gradx_ax
global gradx_sag
global gradx_cor
global grady_ax
global grady_sag
global grady_cor
global axial_m1
global sag_m1
global cor_m1
global scalar

mesh1 = mesh0_X(1:length(mesh0_X)/2);
mesh2 = mesh0_X(length(mesh0_X)/2 + 1 :end );

mesh11 = reshape(mesh1',2,length(mesh0_X)/4);
mesh21 = reshape(mesh2',2,length(mesh0_X)/4);

tri1 = TriRep(source_tri.Triangulation, [mesh11' source_tri.X(:,3)]); 
tri2 = TriRep(source_tri.Triangulation, [source_tri.X(:,1) mesh21']);

% tri1 = TriRep( source_tri.Triangulation,[mesh0_X(1:size( source_tri.X,1),:) source_tri.X(:,3)]); 
% tri2 = TriRep( source_tri.Triangulation,[source_tri.X(:,1) mesh0_X(size( source_tri.X,1)+1:2 * size( source_tri.X,1),:)]); 
% tri3 = TriRep( source_tri.Triangulation,[mesh0_X(2 * size( source_tri.X,1)+1:end,1) source_tri.X(:,2) mesh0_X(2 * size( source_tri.X,1)+1:end,2)]); 

% tri1 = TriRep( source_tri.Triangulation,mesh0_X(1:size( source_tri.X,1),:) ); 
% tri2 = TriRep( source_tri.Triangulation,mesh0_X(size( source_tri.X,1)+1:2 * size( source_tri.X,1),:)); 
% tri3 = TriRep( source_tri.Triangulation,mesh0_X(2 * size( source_tri.X,1)+1:end,:) ); 

% First intersection
b2c_ncoord1_int1 = baryToCart(tri1,  current_tr_int1,  c2b_coord_int1);
b2c_ncoord2_int1 = baryToCart(tri2,  current_tr_int1,  c2b_coord_int1);

lambda = .01;
sqrt_lambda = lambda; % sqrt(lambda);

n1 = size(var_array1,1);

incr = size(source_tri.X,1);
total_var = 4 * incr; % 9
xyz = 2; % 3
% incr_row = in.n + in.size_source;

F1  = zeros(size(var_array1,1),1); 

r_ax  = size(vol_ax, 1);
c_ax  = size(vol_ax, 2);
r_sag = size(vol_sag,1);
c_sag = size(vol_sag,2);


%% First intersection %%
for i = 1: n1
    
    %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1_ax = axial_m1{sub_3_int1(i)} * [b2c_ncoord1_int1(i,:) 1]'; % 3D point to 2D point in the frame coordinates axial_m1_par1{i}
    
    fl  = floor(tmp_v1_ax(1) + 1);
    fl2 = floor(tmp_v1_ax(2) + 1);
    cl  = ceil(tmp_v1_ax(1) + 1);
    cl2 = ceil(tmp_v1_ax(2) + 1);
    
    min_max_r1ax = min(max(fl2,1),r_ax);
    min_max_r2ax = min(max(cl2,1),r_ax);
    min_max_c1ax = min(max(fl,1), c_ax);
    min_max_c2ax = min(max(cl,1), c_ax);
    
    neig = [vol_ax(min_max_r1ax, min_max_c1ax, sub_3_int1(i))   vol_ax(min_max_r1ax, min_max_c2ax, sub_3_int1(i));...
            vol_ax(min_max_r2ax, min_max_c1ax, sub_3_int1(i))   vol_ax(min_max_r2ax, min_max_c2ax, sub_3_int1(i))];
    
    new_im_ax = bilinear_interpolation(tmp_v1_ax(2), tmp_v1_ax(1), double(neig));
    
    %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1_sg = sag_m1{sub_2_int1(i)} * [b2c_ncoord2_int1(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k} sag_m1_par1{i}

    fl  = floor(tmp_v1_sg(1) + 1);
    fl2 = floor(tmp_v1_sg(2) + 1);
    cl  = ceil(tmp_v1_sg(1) + 1);
    cl2 = ceil(tmp_v1_sg(2) + 1);

    min_max_r1sg = min(max(fl2,1),r_sag);
    min_max_r2sg = min(max(cl2,1),r_sag);
    min_max_c1sg = min(max(fl,1), c_sag);
    min_max_c2sg = min(max(cl,1), c_sag);
    
    neig = [vol_sag(min_max_r1sg, min_max_c1sg, sub_2_int1(i))   vol_sag(min_max_r1sg, min_max_c2sg, sub_2_int1(i));...
            vol_sag(min_max_r2sg, min_max_c1sg, sub_2_int1(i))   vol_sag(min_max_r2sg, min_max_c2sg, sub_2_int1(i))];
    
    new_im_sag = bilinear_interpolation(tmp_v1_sg(2), tmp_v1_sg(1), double(neig));
     
    %% Function
    F1(i) =  new_im_ax - new_im_sag;
    
end

F = sum(F1.^2);

if nargout > 1
    
    Js_int1 = zeros(size(var_array1,1), total_var);

    incr_col = xyz * incr;
    
    %% First intersection %%
    for i = 1: n1
        
        %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1_ax = axial_m1{sub_3_int1(i)} * [b2c_ncoord1_int1(i,:) 1]'; % 3D point to 2D point in the frame coordinates axial_m1_par1{i}
        
        fl  = floor(tmp_v1_ax(1) + 1);
        fl2 = floor(tmp_v1_ax(2) + 1);
        cl  = ceil(tmp_v1_ax(1) + 1);
        cl2 = ceil(tmp_v1_ax(2) + 1);
        
        min_max_r1ax = min(max(fl2,1),r_ax);
        min_max_r2ax = min(max(cl2,1),r_ax);
        min_max_c1ax = min(max(fl,1), c_ax);
        min_max_c2ax = min(max(cl,1), c_ax);
        
        %% Gradient in axial view
        
        neigax_gradx = [gradx_ax(min_max_r1ax, min_max_c1ax, sub_3_int1(i)) gradx_ax(min_max_r1ax, min_max_c2ax, sub_3_int1(i));...
                        gradx_ax(min_max_r2ax, min_max_c1ax, sub_3_int1(i)) gradx_ax(min_max_r2ax, min_max_c2ax, sub_3_int1(i))];
        
        neigax_grady = [grady_ax(min_max_r1ax, min_max_c1ax, sub_3_int1(i)) grady_ax(min_max_r1ax, min_max_c2ax, sub_3_int1(i));...
                        grady_ax(min_max_r2ax, min_max_c1ax, sub_3_int1(i)) grady_ax(min_max_r2ax, min_max_c2ax, sub_3_int1(i))];
        
        grad_ax = [-bilinear_interpolation(tmp_v1_ax(2), tmp_v1_ax(1), double(neigax_gradx))  -bilinear_interpolation(tmp_v1_ax(2), tmp_v1_ax(1), double(neigax_grady))];
        
        %     grad_ax  = [gradx_ax( min_max_r, min_max_c, sub_3_int1(i)) grady_ax( min_max_r, min_max_c, sub_3_int1(i))];
        
        %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1_sg = sag_m1{sub_2_int1(i)} * [b2c_ncoord2_int1(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k} sag_m1_par1{i}
        
        fl  = floor(tmp_v1_sg(1) + 1);
        fl2 = floor(tmp_v1_sg(2) + 1);
        cl  = ceil(tmp_v1_sg(1) + 1);
        cl2 = ceil(tmp_v1_sg(2) + 1);
        
        min_max_r1sg = min(max(fl2,1),r_sag);
        min_max_r2sg = min(max(cl2,1),r_sag);
        min_max_c1sg = min(max(fl,1), c_sag);
        min_max_c2sg = min(max(cl,1), c_sag);
        
        %% Gradient in sagittal view
        neigsg_gradx = [gradx_sag(min_max_r1sg, min_max_c1sg, sub_2_int1(i)) gradx_sag(min_max_r1sg, min_max_c2sg, sub_2_int1(i));...
                        gradx_sag(min_max_r2sg, min_max_c1sg, sub_2_int1(i)) gradx_sag(min_max_r2sg, min_max_c2sg, sub_2_int1(i))];
        
        neigsg_grady = [grady_sag(min_max_r1sg, min_max_c1sg, sub_2_int1(i)) grady_sag(min_max_r1sg, min_max_c2sg, sub_2_int1(i));...
                        grady_sag(min_max_r2sg, min_max_c1sg, sub_2_int1(i)) grady_sag(min_max_r2sg, min_max_c2sg, sub_2_int1(i))];
        
        grad_sag = [-bilinear_interpolation(tmp_v1_sg(2), tmp_v1_sg(1), double(neigsg_gradx))  -bilinear_interpolation(tmp_v1_sg(2), tmp_v1_sg(1), double(neigsg_grady))];
        
        %% Jacobian
        w  = zeros(3, total_var);
        w2 = zeros(3, total_var);
        
        ind_w = (tri2.Triangulation(current_tr_int1(i),:)-1)*xyz + 1;
        
        w( 1,ind_w ) = c2b_coord_int1(i,:);
        w2(1,:)      = circshift(w(1,:)', xyz * incr);
        
        w(2,:)  = circshift(w(1,:)', 1);
%         w(3,:)  = circshift(w(1,:)', 2);

        w2(2,:) = w2(1,:);%circshift(w2(1,:)',1);
        w2(3,:) = circshift(w2(1,:)',1); % 2
        w2(1,:) = 0;
            
        Js_int1(i,:) =  (2 * F1(i)) .*( grad_ax * axial_m1{sub_3_int1(i)}(1:2,1:3) * w - grad_sag * sag_m1{sub_2_int1(i)}(1:2,1:3) * w2); %% NOTE:change the indexes of the matrix for 3 components !!!
    end
    
    

    J =  sum(Js_int1);
end




