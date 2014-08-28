function [F J] = myfun4(mesh0_X, in)

mesh1 = mesh0_X(1:length(mesh0_X)/2);
mesh2 = mesh0_X(length(mesh0_X)/2 + 1 :end);

mesh11 = reshape(mesh1',3,length(mesh0_X)/6);
mesh21 = reshape(mesh2',3,length(mesh0_X)/6);

tri1 = TriRep(in.source_tri.Triangulation, mesh11'); % define the new mesh for the axial      mesh0_X(:,1:in.size_source)  mesh0_X(:,1:size(mesh0_X,2)/2)' mesh11'
tri2 = TriRep(in.source_tri.Triangulation, mesh21'); % define the new mesh for the sagittal  mesh21' mesh0_X(:,in.size_source+1:end)  mesh0_X(:,size(mesh0_X,2)/2+1:end)' mesh21'


b2c_ncoord1 = baryToCart(tri1, in.current_tr, in.c2b_coord);
b2c_ncoord2 = baryToCart(tri2, in.current_tr, in.c2b_coord);

incr_row = in.n + in.size_source;
incr_col = 3 * in.size_source;

sub_3 = in.sub_3;
sub_2 = in.sub_2;
rows = in.rows;
cols = in.cols;
vol_ax = in.vol_ax;
vol_sag = in.vol_sag;

F = zeros(in.n, 1);


% F = zeros(in.n + 2 * in.size_source, 1);
% Js = zeros(in.n + 2 * in.size_source, 6 * in.size_source);
%  F = zeros(in.n, 1);
 
source_tri = in.source_tri;
current_tr = in.current_tr;
c2b_coord = in.c2b_coord;
axi_m1_par = in.axial_m1_par;
sag_m1_par = in.sag_m1_par;



% neigax = in.neig_ax;
% neigsg = in.neig_sag;

for i = 1 : in.n
    
    %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    tmp_v1_ax = axi_m1_par{i} * [b2c_ncoord1(i,:) 1]'; % 3D point to 2D point in the frame coordinates axial_m1{sub_3(i)}
    
    %% Make sure the indexes are correct and use bilinear interpolation
    
    fl1ax = floor(tmp_v1_ax(1) + 1);
    fl2ax = floor(tmp_v1_ax(2) + 1);
    cl1ax = ceil(tmp_v1_ax(1) + 1);
    cl2ax = ceil(tmp_v1_ax(2) + 1);
    
    min_max_r1ax = min(max(fl2ax,1),rows);
    min_max_r2ax = min(max(cl2ax,1),rows);
    min_max_c1ax = min(max(fl1ax,1),cols);
    min_max_c2ax = min(max(cl1ax,1),cols);
    
    neigax = [vol_ax(min_max_r1ax, min_max_c1ax, sub_3(i)) vol_ax(min_max_r1ax, min_max_c2ax, sub_3(i));...
              vol_ax(min_max_r2ax, min_max_c1ax, sub_3(i)) vol_ax(min_max_r2ax, min_max_c2ax, sub_3(i))];
    
    new_im_ax = bilinear_interpolation(tmp_v1_ax(2), tmp_v1_ax(1), double(neigax));
    
%     new_im_ax = interp2(in.a_x, in.a_y, vol_ax(:,:,sub_3(i)), tmp_v1_ax(1), tmp_v1_ax(2), 'linear', 0);
    %'ax'    
    %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    tmp_v1_sg = sag_m1_par{i} * [b2c_ncoord2(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k} sag_m1{sub_2(i)}
    
    %% Make sure the indexes are correct
    
    fl1sg = floor(tmp_v1_sg(1) + 1);
    fl2sg = floor(tmp_v1_sg(2) + 1);
    cl1sg = ceil(tmp_v1_sg(1) + 1);
    cl2sg = ceil(tmp_v1_sg(2) + 1);
    
    min_max_r1sg = min(max(fl2sg,1),rows);
    min_max_r2sg = min(max(cl2sg,1),rows);
    min_max_c1sg = min(max(fl1sg,1),cols);
    min_max_c2sg = min(max(cl1sg,1),cols);
    
    neigsg = [vol_sag(min_max_r1sg, min_max_c1sg, sub_2(i)) vol_sag(min_max_r1sg, min_max_c2sg,sub_2(i));...
              vol_sag(min_max_r2sg, min_max_c1sg, sub_2(i)) vol_sag(min_max_r2sg, min_max_c2sg,sub_2(i))];
    
    new_im_sag = bilinear_interpolation(tmp_v1_sg(2), tmp_v1_sg(1), double(neigsg));

%     new_im_sag = interp2(in.a_x, in.a_y, vol_sag(:,:,sub_2(i)), tmp_v1_sg(1), tmp_v1_sg(2), 'linear', 0);
%     'sag'
    %% Function
    
    F(i) = new_im_ax - new_im_sag;

end

F = F./in.n;
% for i = 1:in.size_source
%     
%     mean_tmp1 = mean(tri1.X(in.list_edges{i},:));
%     mean_tmp2 = mean(tri2.X(in.list_edges{i},:));
%     
%     lapl_tri1 = abs(tri1.X(i,:) - mean_tmp1);
%     lapl_tri2 = abs(tri2.X(i,:) - mean_tmp2);
%     
%     norm2_lapl_tri1 = sum( lapl_tri1 );
%     norm2_lapl_tri2 = sum( lapl_tri2 );
%     
%     F(i + in.n)     = in.sqrt_lambda * norm2_lapl_tri1;
%     F(i + incr_row) = in.sqrt_lambda * norm2_lapl_tri2;
% 
% end


if nargout > 1
    
%     J = zeros(in.n + 2*in.size_source, 6 * in.size_source);
    
    for i = 1 : in.n
        
        %% axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1_ax = axi_m1_par{i} * [b2c_ncoord1(i,:) 1]'; % 3D point to 2D point in the frame coordinates axial_m1{sub_3(i)}
        
        %% Make sure the indexes are correct and use bilinear interpolation
        
        fl1ax = floor(tmp_v1_ax(1) + 1);
        fl2ax = floor(tmp_v1_ax(2) + 1);
        cl1ax = ceil(tmp_v1_ax(1) + 1);
        cl2ax = ceil(tmp_v1_ax(2) + 1);
        
        min_max_r1ax = min(max(fl2ax,1),rows);
        min_max_r2ax = min(max(cl2ax,1),rows);
        min_max_c1ax = min(max(fl1ax,1),cols);
        min_max_c2ax = min(max(cl1ax,1),cols);

        neigax_gradx = [in.gradx_ax(min_max_r1ax, min_max_c1ax, sub_3(i)) in.gradx_ax(min_max_r1ax, min_max_c2ax,sub_3(i));...
                        in.gradx_ax(min_max_r2ax, min_max_c1ax, sub_3(i)) in.gradx_ax(min_max_r2ax, min_max_c2ax,sub_3(i))];
        
        neigax_grady = [in.grady_ax(min_max_r1ax, min_max_c1ax, sub_3(i)) in.grady_ax(min_max_r1ax, min_max_c2ax,sub_3(i));...
                        in.grady_ax(min_max_r2ax, min_max_c1ax, sub_3(i)) in.grady_ax(min_max_r2ax, min_max_c2ax,sub_3(i))];
        
        grad_ax = [-bilinear_interpolation(tmp_v1_ax(2), tmp_v1_ax(1), double(neigax_gradx))  -bilinear_interpolation(tmp_v1_ax(2), tmp_v1_ax(1), double(neigax_grady))];
        
        %% sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        tmp_v1_sg = sag_m1_par{i} * [b2c_ncoord2(i,:) 1]'; % 3D point to 2D point in the frame coordinates var_cell{i,j,k} sag_m1{sub_2(i)}
        
        %% Make sure the indexes are correct
        
        fl1sg = floor(tmp_v1_sg(1) + 1);
        fl2sg = floor(tmp_v1_sg(2) + 1);
        cl1sg = ceil(tmp_v1_sg(1) + 1);
        cl2sg = ceil(tmp_v1_sg(2) + 1);
        
        min_max_r1sg = min(max(fl2sg,1),rows);
        min_max_r2sg = min(max(cl2sg,1),rows);
        min_max_c1sg = min(max(fl1sg,1),cols);
        min_max_c2sg = min(max(cl1sg,1),cols);
        
        neigsg_gradx = [in.gradx_sag(min_max_r1sg, min_max_c1sg, sub_2(i)) in.gradx_sag(min_max_r1sg, min_max_c2sg,sub_2(i));...
                        in.gradx_sag(min_max_r2sg, min_max_c1sg, sub_2(i)) in.gradx_sag(min_max_r2sg, min_max_c2sg,sub_2(i))];
        
        neigsg_grady = [in.grady_sag(min_max_r1sg, min_max_c1sg, sub_2(i)) in.grady_sag(min_max_r1sg, min_max_c2sg,sub_2(i));...
                        in.grady_sag(min_max_r2sg, min_max_c1sg, sub_2(i)) in.grady_sag(min_max_r2sg, min_max_c2sg,sub_2(i))];
        
        grad_sag = [-bilinear_interpolation(tmp_v1_sg(2), tmp_v1_sg(1), double(neigsg_gradx))  -bilinear_interpolation(tmp_v1_sg(2), tmp_v1_sg(1), double(neigsg_grady))];

        
        %% Jacobian of warp function qi = w(pi,X), pi is a given data point, and X is the mesh
        
        w  = zeros(3, 2 * incr_col);
        w2 = zeros(3, 2 * incr_col);
        
        %% Get the index of the points in the mesh (tetrahedron) that contains point i
        %% We multiply by 3 because we have 3 coordinates per each point
        ind_w1 = (tri1.Triangulation(current_tr(i),:)-1) * 3 + 1;
        
        w(1,ind_w1 ) = c2b_coord(i,:);
        w2(1,:)      = circshift(w(1,:)', incr_col);
        
        w(2,:)  = circshift(w(1,:)', 1);
        w(3,:)  = circshift(w(1,:)', 2);

        w2(2,:) = circshift(w2(1,:)',1);
        w2(3,:) = circshift(w2(1,:)',2);
        

        %% Jacobian
                
        J(i,:) =  (grad_ax * axi_m1_par{i}(1:2,1:3) * w) - (grad_sag * sag_m1_par{i}(1:2,1:3) * w2);  

    end
    J = J./in.n;
%     for i = 1:in.size_source
%         
%         jsi1 = (i - 1) * 3 + 1;
%         jsi2 = jsi1 + incr_col;
%         J(i + in.n    , jsi1:jsi1+2) = in.sqrt_lambda;%.* lapl_tri1; %(in.sqrt_lambda * ( 1 / norm_lapl_tri1 ) ) .* lapl_tri1;
%         J(i + incr_row, jsi2:jsi2+2) = in.sqrt_lambda;%.* lapl_tri2; %(in.sqrt_lambda * ( 1 / norm_lapl_tri2 ) ) .* lapl_tri2;
%         
%         for j = 1:length(in.list_edges{i})
%             
%             edgein1 = (in.list_edges{i}(j)-1) * 3 + 1;
%             edgein2 = edgein1 + incr_col;
%             J(i + in.n,     edgein1:edgein1+2) = J(i + in.n,     jsi1:jsi1+2).*[in.scalar(i) in.scalar(i) in.scalar(i)];
%             J(i + incr_row, edgein2:edgein2+2) = J(i + incr_row, jsi2:jsi2+2).*[in.scalar(i) in.scalar(i) in.scalar(i)];
%             
%         end
%         
%     end

end




