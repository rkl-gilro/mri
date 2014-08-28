function [F, J] = myfun_unc_orthoJp(mesh0_X)


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

mesh1 = mesh0_X(1:length(mesh0_X)/3);
mesh2 = mesh0_X(length(mesh0_X)/3 + 1 :2 * length(mesh0_X)/3 );
mesh3 = mesh0_X(2*length(mesh0_X)/3 + 1 :end );

mesh11 = reshape(mesh1',2,length(mesh0_X)/6);
mesh21 = reshape(mesh2',2,length(mesh0_X)/6);
mesh31 = reshape(mesh3',2,length(mesh0_X)/6);


tri1 = TriRep(source_tri.Triangulation, [mesh11' source_tri.X(:,3)]); 
tri2 = TriRep(source_tri.Triangulation, [source_tri.X(:,1) mesh21']);
tri3 = TriRep(source_tri.Triangulation, [mesh31(1,:)' source_tri.X(:,2) mesh31(2,:)']);

% tri1 = TriRep( source_tri.Triangulation,[mesh0_X(1:size( source_tri.X,1),:) source_tri.X(:,3)]); 
% tri2 = TriRep( source_tri.Triangulation,[source_tri.X(:,1) mesh0_X(size( source_tri.X,1)+1:2 * size( source_tri.X,1),:)]); 
% tri3 = TriRep( source_tri.Triangulation,[mesh0_X(2 * size( source_tri.X,1)+1:end,1) source_tri.X(:,2) mesh0_X(2 * size( source_tri.X,1)+1:end,2)]); 

% tri1 = TriRep( source_tri.Triangulation,mesh0_X(1:size( source_tri.X,1),:) ); 
% tri2 = TriRep( source_tri.Triangulation,mesh0_X(size( source_tri.X,1)+1:2 * size( source_tri.X,1),:)); 
% tri3 = TriRep( source_tri.Triangulation,mesh0_X(2 * size( source_tri.X,1)+1:end,:) ); 

% First intersection
b2c_ncoord1_int1 = baryToCart(tri1,  current_tr_int1,  c2b_coord_int1);
b2c_ncoord2_int1 = baryToCart(tri2,  current_tr_int1,  c2b_coord_int1);

% Second intersection
b2c_ncoord1_int2 = baryToCart(tri1,  current_tr_int2,  c2b_coord_int2);
b2c_ncoord2_int2 = baryToCart(tri3,  current_tr_int2,  c2b_coord_int2);

% Third intersection
b2c_ncoord1_int3 = baryToCart(tri2,  current_tr_int3,  c2b_coord_int3);
b2c_ncoord2_int3 = baryToCart(tri3,  current_tr_int3,  c2b_coord_int3);

lambda = .01;
sqrt_lambda = sqrt(lambda); % lambda; % 

n1 = size(var_array1,1);
n2 = size(var_array2,1);
n3 = size(var_array3,1);

incr = size(source_tri.X,1);
total_var = 6 * incr; % 9
xyz = 2; % 3
% incr_row = in.n + in.size_source;

% F1  = zeros(size(var_array1,1),1); % + size(source_tri.X,1) + size(source_tri.X,1),1);
% F2  = zeros(size(var_array2,1),1); 
% F3  = zeros(size(var_array3,1),1); 
% F_s = zeros(size(mesh0_X,1),1);

r_ax  = size(vol_ax, 1);
c_ax  = size(vol_ax, 2);
r_sag = size(vol_sag,1);
c_sag = size(vol_sag,2);
r_cor = size(vol_cor,1);
c_cor = size(vol_cor,2);

funcs = {@Intersection1, @Intersection2, @Intersection3};
arguments = {n1, b2c_ncoord1_int1, axial_m1, sub_3_int1,  vol_ax, b2c_ncoord2_int1, sag_m1, sub_2_int1, vol_sag;...
             n2, b2c_ncoord1_int2, axial_m1, sub_3_int2,  vol_ax, b2c_ncoord2_int2, cor_m1, sub_2_int2, vol_cor;...
             n3, b2c_ncoord1_int3,   sag_m1, sub_2_int3, vol_sag, b2c_ncoord2_int3, cor_m1, sub_3_int3, vol_cor};
solutions = cell(1,3);

parfor i = 1:3
    solutions{i} = funcs{i}(arguments{i,:});
end
F1 = solutions{1};
F2 = solutions{2};
F3 = solutions{3};

%% Smooth term %%
for i = 1:size( source_tri.X,1)
    
    mean_tmp1 = mean(tri1.X(list_edges{i},:));
    mean_tmp2 = mean(tri2.X(list_edges{i},:));
    mean_tmp3 = mean(tri3.X(list_edges{i},:));
    
    lapl_tri1 = abs(tri1.X(i,:) - mean_tmp1);
    lapl_tri2 = abs(tri2.X(i,:) - mean_tmp2);
    lapl_tri3 = abs(tri3.X(i,:) - mean_tmp3);
    
    %% Functional %%
    F_s( i )         =  sqrt_lambda * sum( lapl_tri1(1:2).^2 ); % sqrt_lambda * (sum( lapl_tri1.^2 ));
    F_s(i + incr)    =  sqrt_lambda * sum( lapl_tri2(2:3).^2 ); % sqrt_lambda * (sum( lapl_tri2.^2 ));
    F_s(i + 2*incr)  =  sqrt_lambda * sum( [lapl_tri3(1) lapl_tri3(3)].^2 ); % sqrt_lambda * (sum( lapl_tri3.^2 ));
    
end
    
F = (.5/(n1 + n2 + n3 + 3*size( source_tri.X,1)) ).* (sum(F1.^2) + sum(F2.^2) + sum(F3.^2) + sum(F_s));


if nargout > 1
    
    Js_int1 = zeros(size(var_array1,1), total_var);
    Js_int2 = zeros(size(var_array2,1), total_var);
    Js_int3 = zeros(size(var_array3,1), total_var);
    Js_s    = zeros(            3*incr, total_var);

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
            
        Js_int1(i,:) = (F1(i)) .*  (grad_ax * axial_m1{sub_3_int1(i)}(1:2,1:3) * w - grad_sag * sag_m1{sub_2_int1(i)}(1:2,1:3) * w2); %% NOTE:change the indexes of the matrix for 3 components !!!
    end
    
    %% Second intersection %%
    for i = 1: n2
        
        %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1_ax = axial_m1{sub_3_int2(i)} * [b2c_ncoord1_int2(i,:) 1]'; % 3D point to 2D point in the frame coordinates axial_m1_par2{i}
        
        fl  = floor(tmp_v1_ax(1) + 1);
        fl2 = floor(tmp_v1_ax(2) + 1);
        cl  = ceil(tmp_v1_ax(1) + 1);
        cl2 = ceil(tmp_v1_ax(2) + 1);
        
        min_max_r1ax = min(max(fl2,1),r_ax);
        min_max_r2ax = min(max(cl2,1),r_ax);
        min_max_c1ax = min(max(fl,1), c_ax);
        min_max_c2ax = min(max(cl,1), c_ax);
        
        %% Gradient in axial view
        
        neigax_gradx = [gradx_ax(min_max_r1ax, min_max_c1ax, sub_3_int2(i)) gradx_ax(min_max_r1ax, min_max_c2ax, sub_3_int2(i));...
                        gradx_ax(min_max_r2ax, min_max_c1ax, sub_3_int2(i)) gradx_ax(min_max_r2ax, min_max_c2ax, sub_3_int2(i))];
        
        neigax_grady = [grady_ax(min_max_r1ax, min_max_c1ax, sub_3_int2(i)) grady_ax(min_max_r1ax, min_max_c2ax, sub_3_int2(i));...
                        grady_ax(min_max_r2ax, min_max_c1ax, sub_3_int2(i)) grady_ax(min_max_r2ax, min_max_c2ax, sub_3_int2(i))];
        
        grad_ax = [-bilinear_interpolation(tmp_v1_ax(2), tmp_v1_ax(1), double(neigax_gradx))  -bilinear_interpolation(tmp_v1_ax(2), tmp_v1_ax(1), double(neigax_grady))];
        
        %% coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1_cr =  cor_m1{sub_2_int2(i)} * [b2c_ncoord2_int2(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}
        
        fl  = floor(tmp_v1_cr(1) + 1);
        fl2 = floor(tmp_v1_cr(2) + 1);
        cl  = ceil(tmp_v1_cr(1) + 1);
        cl2 = ceil(tmp_v1_cr(2) + 1);
        
        min_max_r1cr = min(max(fl2,1),r_cor);
        min_max_r2cr = min(max(cl2,1),r_cor);
        min_max_c1cr = min(max(fl,1), c_cor);
        min_max_c2cr = min(max(cl,1), c_cor);
        
        %% Gradient in coronal view
        
        neigcr_gradx = [gradx_cor(min_max_r1cr, min_max_c1cr, sub_2_int2(i)) gradx_cor(min_max_r1cr, min_max_c2cr, sub_2_int2(i));...
                        gradx_cor(min_max_r2cr, min_max_c1cr, sub_2_int2(i)) gradx_cor(min_max_r2cr, min_max_c2cr, sub_2_int2(i))];
        
        neigcr_grady = [grady_cor(min_max_r1cr, min_max_c1cr, sub_2_int2(i)) grady_cor(min_max_r1cr, min_max_c2cr, sub_2_int2(i));...
                        grady_cor(min_max_r2cr, min_max_c1cr, sub_2_int2(i)) grady_cor(min_max_r2cr, min_max_c2cr, sub_2_int2(i))];
        
        grad_cor = [-bilinear_interpolation(tmp_v1_cr(2), tmp_v1_cr(1), double(neigcr_gradx))  -bilinear_interpolation(tmp_v1_cr(2), tmp_v1_cr(1), double(neigcr_grady))];
                
        %% Jacobian
        w  = zeros(3, total_var);
        w2 = zeros(3, total_var);
        
        ind_w = (tri1.Triangulation(current_tr_int2(i),:)-1) * xyz + 1;
        
        w( 1,ind_w ) = c2b_coord_int2(i,:);
        w2(1,:)      = circshift(w(1,:)', 2 * xyz * incr );
        
        w( 2,:) = circshift(w(1,:)', 1);
%         w( 3,:) = circshift(w(1,:)', 2);
%         w2(2,:) = circshift(w2(1,:)',1);
        w2(3,:) = circshift(w2(1,:)',1); % 2
        
        Js_int2(i,:) = (F2(i)) .*  ((grad_ax * axial_m1{sub_3_int2(i)}(1:2,1:3) * w) - (grad_cor * cor_m1{sub_2_int2(i)}(1:2,1:3) * w2));
        
    end
    %% Third intersection %%
    for i = 1: n3
        
        %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1_sg = sag_m1{sub_2_int3(i)} * [b2c_ncoord1_int3(i,:) 1]'; % 3D point to 2D point in the frame coordinates sag_m1_par2{i}
        
        fl  = floor(tmp_v1_sg(1) + 1);
        fl2 = floor(tmp_v1_sg(2) + 1);
        cl  = ceil(tmp_v1_sg(1)  + 1);
        cl2 = ceil(tmp_v1_sg(2)  + 1);
        
        min_max_r1sg = min(max(fl2,1),r_sag);
        min_max_r2sg = min(max(cl2,1),r_sag);
        min_max_c1sg = min(max(fl,1), c_sag);
        min_max_c2sg = min(max(cl,1), c_sag);
        
        %% Gradient in sagittal view
        
        neigsg_gradx = [gradx_sag(min_max_r1sg, min_max_c1sg, sub_2_int3(i)) gradx_sag(min_max_r1sg, min_max_c2sg, sub_2_int3(i));...
                        gradx_sag(min_max_r2sg, min_max_c1sg, sub_2_int3(i)) gradx_sag(min_max_r2sg, min_max_c2sg, sub_2_int3(i))];
        
        neigsg_grady = [grady_sag(min_max_r1sg, min_max_c1sg, sub_2_int3(i)) grady_sag(min_max_r1sg, min_max_c2sg, sub_2_int3(i));...
                        grady_sag(min_max_r2sg, min_max_c1sg, sub_2_int3(i)) grady_sag(min_max_r2sg, min_max_c2sg, sub_2_int3(i))];
        
        grad_sag = [-bilinear_interpolation(tmp_v1_sg(2), tmp_v1_sg(1), double(neigsg_gradx))  -bilinear_interpolation(tmp_v1_sg(2), tmp_v1_sg(1), double(neigsg_grady))];
        
        %% coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1_cr = cor_m1{sub_3_int3(i)} * [b2c_ncoord2_int3(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k}cor_m1_par2{i}
        
        fl  = floor(tmp_v1_cr(1) + 1);
        fl2 = floor(tmp_v1_cr(2) + 1);
        cl  = ceil(tmp_v1_cr(1)  + 1);
        cl2 = ceil(tmp_v1_cr(2)  + 1);
        
        min_max_r1cr = min(max(fl2,1),r_cor);
        min_max_r2cr = min(max(cl2,1),r_cor);
        min_max_c1cr = min(max(fl,1), c_cor);
        min_max_c2cr = min(max(cl,1), c_cor);
        
        %% Gradient in coronal view
        
        neigcr_gradx = [gradx_cor(min_max_r1cr, min_max_c1cr, sub_3_int3(i)) gradx_cor(min_max_r1cr, min_max_c2cr, sub_3_int3(i));...
                        gradx_cor(min_max_r2cr, min_max_c1cr, sub_3_int3(i)) gradx_cor(min_max_r2cr, min_max_c2cr, sub_3_int3(i))];
        
        neigcr_grady = [grady_cor(min_max_r1cr, min_max_c1cr, sub_3_int3(i)) grady_cor(min_max_r1cr, min_max_c2cr, sub_3_int3(i));...
                        grady_cor(min_max_r2cr, min_max_c1cr, sub_3_int3(i)) grady_cor(min_max_r2cr, min_max_c2cr, sub_3_int3(i))];
        
        grad_cor = [-bilinear_interpolation(tmp_v1_cr(2), tmp_v1_cr(1), double(neigcr_gradx))  -bilinear_interpolation(tmp_v1_cr(2), tmp_v1_cr(1), double(neigcr_grady))];

        %% Jacobian
        w_tmp  = zeros(3, total_var);
        w  = zeros(3, total_var);
        w2 = zeros(3, total_var);
        
        ind_w = (tri3.Triangulation(current_tr_int3(i),:)-1) * xyz + 1;
        
        w_tmp( 1,ind_w ) = c2b_coord_int3(i,:);
        w(1,:)   = circshift(w_tmp(1,:)', xyz * incr);
        w2(1,:)  = circshift(w_tmp(1,:)', 2 * xyz * incr);
        
        w(2,:) = w(1,:);%circshift(w(1,:)',  1);
        w(3,:) = circshift(w(1,:)',  1); % 2
        w(1,:) = 0;
        
%         w2(2,:) = circshift(w2(1,:)',1);
        w2(3,:) = circshift(w2(1,:)',1); % 2
        
        Js_int3(i,:) = (F3(i)) .* ( (grad_cor * cor_m1{sub_3_int3(i)}(1:2,1:3) * w2) - (grad_sag * sag_m1{sub_2_int3(i)}(1:2,1:3) * w)); %
        
    end
    
    %% Smooth term %%
    for i = 1:size( source_tri.X,1)
        
        mean_tmp1 = mean(tri1.X(list_edges{i},:));
        mean_tmp2 = mean(tri2.X(list_edges{i},:));
        mean_tmp3 = mean(tri3.X(list_edges{i},:));
        
%         lapl_tri1 = abs(tri1.X(i,:) - mean_tmp1);
%         lapl_tri2 = abs(tri2.X(i,:) - mean_tmp2);
%         lapl_tri3 = abs(tri3.X(i,:) - mean_tmp3);

        lapl_tri1 = (tri1.X(i,:) - mean_tmp1);
        lapl_tri2 = (tri2.X(i,:) - mean_tmp2);
        lapl_tri3 = (tri3.X(i,:) - mean_tmp3);
        
        %% Gradient %%
        jsi1 = (i - 1) * xyz + 1;
        jsi2 = jsi1 + incr_col;
        jsi3 = jsi2 + incr_col;
        
        Js_s(i    ,        jsi1:jsi1+1) =  (lambda).* lapl_tri1(1:2); % 
        Js_s(i + incr,     jsi2:jsi2+1) =  (lambda).* lapl_tri2(2:3); % sqrt_lambda; % 
        Js_s(i + 2 * incr, jsi3:jsi3+1) =  (lambda).* [lapl_tri3(1) lapl_tri3(3)]; % xsqrt_lambda; % 
       
        
        for j = 1:length(list_edges{i})
            
            edgein1 = (list_edges{i}(j)-1) * xyz + 1;
            edgein2 = edgein1 + incr_col;
            edgein3 = edgein2 + incr_col;
            
            Js_s(i   ,         edgein1:edgein1+1) = Js_s(i ,         jsi1:jsi1+1).*[scalar(i) scalar(i)];
            Js_s(i + incr,     edgein2:edgein2+1) = Js_s(i + incr,   jsi2:jsi2+1).*[scalar(i) scalar(i)];
            Js_s(i + 2 * incr, edgein3:edgein3+1) = Js_s(i + 2*incr, jsi3:jsi3+1).*[scalar(i) scalar(i)];
            
%             Js_s(i   ,         edgein1:edgein1+2) = Js_s(i ,         jsi1:jsi1+2).*[scalar(i) scalar(i) scalar(i)];
%             Js_s(i + incr,     edgein2:edgein2+2) = Js_s(i + incr,   jsi2:jsi2+2).*[scalar(i) scalar(i) scalar(i)];
%             Js_s(i + 2 * incr, edgein3:edgein3+2) = Js_s(i + 2*incr, jsi3:jsi3+2).*[scalar(i) scalar(i) scalar(i)];
            
        end
        
    end

    J = (1/(n1 + n2 + n3 + 3*size( source_tri.X,1)) ).*(sum(Js_int1) + sum(Js_int2) + sum(Js_int3) + sum(Js_s));
end




