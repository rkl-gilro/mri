


%% Calculate the segmented contours in RCS, the originals and after registration (new_ & cont)

%% Axial new points from the segmentation
new_ax_p = cell(size(vol_ax,3),1);
ax_cont  = cell(size(vol_ax,3),1);

for slices = 1:size(vol_ax,3)
    [tmp1, tmp2] = new_contours(vol_ax, axial_m, axial_m1, source_tri, target_tri_ax,  vert_axial, slices);
    new_ax_p{slices} = tmp1;
    ax_cont{slices}  = tmp2;
end

%% Sagittal new points from the segmentation
new_sag_p = cell(size(vol_sag,3),1);
sag_cont  = cell(size(vol_sag,3),1);

for slices = 1:size(vol_sag,3)
    [tmp1, tmp2] = new_contours(vol_sag, sag_m, sag_m1, source_tri, target_tri_sag,  vert_sagittal, slices);
    new_sag_p{slices} = tmp1;
    sag_cont{slices}  = tmp2;
end

%% Coronal new points from the segmentation
new_cor_p = cell(size(vol_cor,3),1);
cor_cont  = cell(size(vol_cor,3),1);

for slices = 1:size(vol_cor,3)
    [tmp1, tmp2] = new_contours(vol_cor, cor_m, cor_m1, source_tri, target_tri_cor,  vert_coronal, slices);
    new_cor_p{slices} = tmp1;
    cor_cont{slices}  = tmp2;    
end

%% Calculate the intersection between planes and the points in the curve contours that belong to it
%% line ax + b

%% First intersection
int1_ax = [9 10];
int1_sag = [9 10];

[X_ax, Y_ax] = meshgrid(int1_ax, int1_sag);

for i=1:length(X_ax(:))
    line_int1(i).a = var_cell1{X_ax(i), Y_ax(i),end};% - var_cell1{X(i), Y(i),1};
    line_int1(i).b = var_cell1{X_ax(i), Y_ax(i),1};
end

%% Second intersection
int2_ax = [9 10];
int2_cor = [9 10];

[X_sag, Y_sag] = meshgrid(int2_ax, int2_cor);

for i=1:length(X_sag(:))
    line_int2(i).a = var_cell2{X_sag(i), Y_sag(i),end};% - var_cell2{X(i), Y(i),1};
    line_int2(i).b = var_cell2{X_sag(i), Y_sag(i),1};
end

%% Third intersection
int3_cor = [9 10];
int3_sag = [9 10];

[X_cor, Y_cor] = meshgrid(int3_cor, int3_sag);

for i=1:length(X_cor(:))
    line_int3(i).a = var_cell3{X_cor(i), Y_cor(i),end};% - var_cell3{X(i), Y(i),1};
    line_int3(i).b = var_cell3{X_cor(i), Y_cor(i),1};
end

%% Calculate the intersection points
% (x - x1) / (x2 - x1) = (y - y1) / (y2 - y1) = (z - z1) / (z2 - z1)
% d = abs(cross(P-Q1,P-Q2))/abs(Q2-Q1);

%% First intersection axial
for i=1:length(int1_ax)

    
    
    for j=1:size(new_ax_p{int1_ax(i)},2)
        if ~isempty(new_ax_p{int1_ax(i)}(:,j))
            d_new_int1_ax{i,j} = norm((cross(new_ax_p{int1_ax(i)}(1:3,j)' - line_int1(i).a,new_ax_p{int1_ax(i)}(1:3,j)' - line_int1(i).b)))/norm((line_int1(i).a - line_int1(i).b));
            d_old_int1_ax{i,j} = norm((cross(    ax_cont{int1_ax(i)}(j,:) - line_int1(i).a,    ax_cont{int1_ax(i)}(j,:) - line_int1(i).b)))/norm((line_int1(i).a - line_int1(i).b));
        end
        
    end

    
end

%% First intersection sagittal
for i=1:length(int1_sag)

    
    if ~isempty(new_sag_p{int1_sag(i)})
        for j=1:size(new_sag_p{int1_sag(i)},2)
            
            d_new_int1_sag{i,j} = norm(cross(new_sag_p{int1_sag(i)}(1:3,j)' - line_int1(i).a,new_sag_p{int1_sag(i)}(1:3,j)' - line_int1(i).b))/norm(line_int1(i).a - line_int1(i).b);
            d_old_int1_sag{i,j} = norm(cross(    sag_cont{int1_sag(i)}(j,:) - line_int1(i).a,    sag_cont{int1_sag(i)}(j,:) - line_int1(i).b))/norm(line_int1(i).a - line_int1(i).b);
            
        end
    end

    
end

%% Compare the intersections before and after registration
for i=1%:length(int1_ax)
    
    %% New
    [~, ind1] = sort([d_new_int1_ax{i,:}]);
    [~, ind2] = sort([d_new_int1_sag{i,:}]);

    norm(new_ax_p{int1_ax(i)}(1:3,ind1(2)) - new_sag_p{int1_ax(i)}(1:3,ind2(1)))
    norm(new_ax_p{int1_ax(i)}(1:3,ind1(1)) - new_sag_p{int1_ax(i)}(1:3,ind2(2)))
    
    %% Old
    [~, ind1] = sort([d_old_int1_ax{i,:}]);
    [~, ind2] = sort([d_old_int1_sag{i,:}]);
    
    norm(ax_cont{int1_ax(i)}(ind1(2),:) - sag_cont{int1_ax(i)}(ind2(2),:))
    norm(ax_cont{int1_ax(i)}(ind1(1),:) - sag_cont{int1_ax(i)}(ind2(1),:))

%     for j=1:size(new_ax_p{int1_ax(i)},2)
%         if ~isempty(new_ax_p{int1_ax(i)}(:,j))
%             d_new_int1_ax{i,j} = sum(abs(cross(line_int1(i).a - line_int1(i).b,new_ax_p{int1_ax(i)}(1:3,j)' - line_int1(i).b)))/sum(abs(line_int1(i).a - line_int1(i).b));
%             d_old_int1_ax{i,j} = sum(abs(cross(line_int1(i).a - line_int1(i).b, ax_cont{int1_ax(i)}(j,:) - line_int1(i).b)))/sum(abs(line_int1(i).a - line_int1(i).b));
%         end
%         
%     end

    
end

% %% Third intersection 
% for i=1:length(int1_ax)
% 
%     
%     if ~isempty(new_ax_p{int1_ax(i)})
%         for j=1:size(new_ax_p{int1_ax(i)},2)
%             
%             d_new{i,j} = sum(abs(cross(line_int1(i).a - line_int1(i).b,new_ax_p{int1_ax(i)}(1:3,j)' - line_int1(i).b)))/sum(abs(line_int1(i).a - line_int1(i).b));
%             d_old{i,j} = sum(abs(cross(line_int1(i).a - line_int1(i).b, ax_cont{int1_ax(i)}(1:3,j)' - line_int1(i).b)))/sum(abs(line_int1(i).a - line_int1(i).b));
%             
%         end
%     end
% 
%     
% end

