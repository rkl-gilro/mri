function [cont_reg, cont] = new_contours(vol, M, M1, source, target,  vert, slices)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Calculate the new contours(registered) given the ones segmented 
%% by the expert(not registered)
%%
%% Inputs:  1. vol -> MxNxS original volume from one view (T2-W)
%%          2. M   -> cell(S,1) containing the transformation from pixel to RCS
%%          3. M1   -> cell(S,1) containing the transformation from RCS to pixel
%%          4. source   -> source triangulation (TriRep) before registration
%%          5. target   -> target triangulation (TriRep) after registration
%%          6. vert   -> cell(S,1) of arrays Tx2 with (i,j) pixel coordinates
%%          7. slices   -> scalar refering to the current slice 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% number of slices
[~,~,N] = size(vol);
cont = [];
cont_reg = [];

% cont = cell(N,1);
% new_cont = cell(N,1);
% tmp_v1 = cell(N,1);
% new_cont2 = cell(N,1);
% tmp_v12 = cell(N,1);

i = slices; %1:N

if ~isempty(vert{i})%.contour{1}
    
    x = [];
    y = [];
    z = [];
    
    for k=1:size(vert{i},1)%.contour{1}
        
        %% NOTE: the contours are given by (j, i) = (x,y) this is why it is written in this order
        p = M{i} * [vert{i}(k,1)-1 vert{i}(k,2)-1 1]';%.contour{1}
        
        x = [x p(1)];
        y = [y p(2)];
        z = [z p(3)];
        
    end
    cont = [x' y' z'];
    
    %% Using the new triangulation calculate the new coordinates
    current_tr = tsearchn(source.X, source.Triangulation,[cont(:,1) cont(:,2) cont(:,3)]); % calculate the tetrahedron where p_3d belongs
    
    c2b_coord  = cartToBary(source, current_tr,[cont(:,1) cont(:,2) cont(:,3)]); % get the barycentric coordinates
    
    new_cont{i}  = baryToCart(target, current_tr, c2b_coord);
    new_cont2{i} = baryToCart(source, current_tr, c2b_coord);
    
    tmp_v1{i}  = M1{i} * [new_cont{i}  ones(size(new_cont{i},1),1)]';
    tmp_v12{i} = M1{i} * [new_cont2{i} ones(size(new_cont2{i},1),1)]';
    
    cont_reg = M{slices} * [tmp_v12{slices}(1,:) - (tmp_v1{slices}(1,:)-tmp_v12{slices}(1,:))-1; tmp_v12{slices}(2,:)- (tmp_v1{slices}(2,:)-tmp_v12{slices}(2,:))-1; ones(1,size(tmp_v1{slices},2))];
else
    cont_reg = [];
    cont = [];
end




