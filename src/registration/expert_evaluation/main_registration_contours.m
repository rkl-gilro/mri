%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Main program for computing the segmented contours from the expert,
%% in the new registered images.
%% If you want to plot the results set 'show' to 1, it is set to 0 (default)
%%
%% First vert_axial, vert_sagittal & vert_coronal have to be defined
%% E.g. for the axial view: 
%% 1. load('/home/raquel/Documents/repositories/mri-segmentation/resources/manual_seg/LASSALAS/mat/uterus/lassalas_uterus_ax.mat')
%% 2. vert_axial = vertex;
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
show = 1;

az = 45; % [45 135 225 -45]
el = -45; % [45 -45]
slices = [9 9 12];

%% Axial new points from the segmentation
[new_ax_p, ax_cont] = new_contours(vol_ax, axial_m, axial_m1, source_tri, target_tri_ax,  vert_axial, slices(1));


%% Sagittal new points from the segmentation
[new_sag_p, sag_cont] = new_contours(vol_sag, sag_m, sag_m1, source_tri, target_tri_sag,  vert_sagittal, slices(2));


%% Coronal new points from the segmentation
[new_cor_p, cor_cont] = new_contours(vol_cor, cor_m, cor_m1, source_tri, target_tri_cor,  vert_coronal, slices(3));

%% Plot of the results, volume with the contours
[rows_ax, cols_ax,~]   = size(vol_ax);
[rows_sag, cols_sag,~] = size(vol_sag);
[rows_cor, cols_cor,~] = size(vol_cor);

if show
    figure;
    % Deformation Axial and Sagittal
    aa(:,:,1) = convert2u8(new_axial(:,:,slices(1)));
    aa(:,:,2) = convert2u8(new_axial(:,:,slices(1)));
    aa(:,:,3) = convert2u8(new_axial(:,:,slices(1)));
    
    [X, Y, Z, triTexture] = compute_triTexture(axial_m{slices(1)},aa, [0 rows_ax-1],[0 cols_ax-1]);
    surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
    hold on
    axis equal
    
    ss(:,:,2) = convert2u8(new_sagittal(:,:,slices(2)));
    ss(:,:,1) = convert2u8(new_sagittal(:,:,slices(2)));
    ss(:,:,3) = convert2u8(new_sagittal(:,:,slices(2)));
    
    [X, Y, Z, triTexture] = compute_triTexture(sag_m{slices(2)},ss, [0 rows_sag-1],[0 cols_sag-1]);
    surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
    hold on
    
    ss(:,:,2) = convert2u8(new_coronal(:,:,slices(3)));
    ss(:,:,1) = convert2u8(new_coronal(:,:,slices(3)));
    ss(:,:,3) = convert2u8(new_coronal(:,:,slices(3)));
    
    [X, Y, Z, triTexture] = compute_triTexture(cor_m{slices(3)},ss, [0 rows_cor-1],[0 cols_cor-1]);
    surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
    hold on
    
    %% Plot the new contours
    % new_ax_p = axial_m{slices(1)} * [tmp_v1_ax2{slices(1)}(1,:) - (tmp_v1_ax{slices(1)}(1,:)-tmp_v1_ax2{slices(1)}(1,:)) - 1; tmp_v1_ax2{slices(1)}(2,:)- (tmp_v1_ax{slices(1)}(2,:)-tmp_v1_ax2{slices(1)}(2,:))- 1; ones(1,size(tmp_v1_ax{slices(1)},2))];
    plot3(new_ax_p(1,11), new_ax_p(2,11), new_ax_p(3,11), 'm*');hold on%, 'Linewidth',3);hold on
    plot3(new_ax_p(1,16), new_ax_p(2,16), new_ax_p(3,16), 'm*');hold on
    % new_sag_p = sag_m{slices(2)} * [tmp_v1_sag2{slices(2)}(1,:) - (tmp_v1_sag{slices(2)}(1,:)-tmp_v1_sag2{slices(2)}(1,:))-1; tmp_v1_sag2{slices(2)}(2,:) - (tmp_v1_sag{slices(2)}(2,:)-tmp_v1_sag2{slices(2)}(2,:))-1; ones(1,size(tmp_v1_sag{slices(2)},2))];
    plot3(new_sag_p(1,11), new_sag_p(2,11), new_sag_p(3,11), 'g*');hold on%, 'Linewidth',3);hold on
    plot3(new_sag_p(1,24), new_sag_p(2,24), new_sag_p(3,24), 'g*');hold on
    % new_cor_p = cor_m{slices(3)} * [tmp_v1_cor2{slices(3)}(1,:) - (tmp_v1_cor{slices(3)}(1,:)-tmp_v1_cor2{slices(3)}(1,:))-1; tmp_v1_cor2{slices(3)}(2,:)- (tmp_v1_cor{slices(3)}(2,:)-tmp_v1_cor2{slices(3)}(2,:))-1; ones(1,size(tmp_v1_cor{slices(3)},2))];
    plot3(new_cor_p(1,:), new_cor_p(2,:), new_cor_p(3,:), 'b', 'Linewidth',3);hold on
    
    title('Registered');
    
    view([az el])
    
    axis equal
    axis off; grid off
    
    [new_ax_p(1,16), new_ax_p(2,16), new_ax_p(3,16)]
    [new_sag_p(1,24), new_sag_p(2,24), new_sag_p(3,24)]
    %% Visualize Original %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    % Original Axial and Sagittal
    
    figure;
    
    aa(:,:,1) = convert2u8(vol_ax(:,:,slices(1)));
    aa(:,:,2) = convert2u8(vol_ax(:,:,slices(1)));
    aa(:,:,3) = convert2u8(vol_ax(:,:,slices(1)));
    
    [X, Y, Z, triTexture] = compute_triTexture(axial_m{slices(1)},aa,[0 rows_ax-1],[0 cols_ax-1]);
    surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
    hold on
    axis equal
    
    ss(:,:,2) = convert2u8(vol_sag(:,:,slices(2)));
    ss(:,:,1) = convert2u8(vol_sag(:,:,slices(2)));
    ss(:,:,3) = convert2u8(vol_sag(:,:,slices(2)));
    
    [X, Y, Z, triTexture] = compute_triTexture(sag_m{slices(2)},ss,[0 rows_sag-1],[0 cols_sag-1]);
    surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
    hold on
    
    ss(:,:,2) = convert2u8(vol_cor(:,:,slices(3)));
    ss(:,:,1) = convert2u8(vol_cor(:,:,slices(3)));
    ss(:,:,3) = convert2u8(vol_cor(:,:,slices(3)));
    
    [X, Y, Z, triTexture] = compute_triTexture(cor_m{slices(3)},ss, [0 rows_cor-1],[0 cols_cor-1]);
    surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
    hold on
    
    %% Plot the original contours
%     plot3( ax_cont(:,1),  ax_cont(:,2),  ax_cont(:,3),'m', 'Linewidth',3);hold on % {slices(1)}
    plot3( ax_cont(12,1),  ax_cont(12,2),  ax_cont(12,3),'m*');hold on
    plot3( ax_cont(16,1),  ax_cont(16,2),  ax_cont(16,3),'m*');hold on
    
%     plot3(sag_cont(:,1), sag_cont(:,2), sag_cont(:,3),'g', 'Linewidth',3);hold on
    plot3(sag_cont(11,1), sag_cont(11,2), sag_cont(11,3),'g*');hold on
    plot3(sag_cont(25,1), sag_cont(25,2), sag_cont(25,3),'g*');hold on
    
    plot3(cor_cont(:,1), cor_cont(:,2), cor_cont(:,3),'b', 'Linewidth',3);hold on
    
    view([az el])
    
    title('Not registered');
    axis equal
    axis off; grid off

end
% clear except- source_tri target_tri_ax target_tri_sag target_tri_cor ...
%                         axial_m axial_m1 sag_m sag_m1 cor_m cor_m1 ... %'new_axial', 'new_sagittal', 'new_coronal',..
%                         vol_ax vol_sag vol_cor views dcmdir ...
%                         optimizer iter measure_errors