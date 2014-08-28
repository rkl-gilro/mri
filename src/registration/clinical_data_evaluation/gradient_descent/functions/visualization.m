function visualization( views, slices )
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Visualisation of the registration results with the original images,
%%  and using the view3d GUI
%%
%%  
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

global vol_ax
global vol_sag
global vol_cor

global new_axial
global new_sagittal
global new_coronal

global axial_m
global sag_m
global cor_m

rows_ax = size(vol_ax,1);
cols_ax = size(vol_ax,2);

rows_sag = size(vol_sag,1);
cols_sag = size(vol_sag,2);

rows_cor = size(vol_cor,1);
cols_cor = size(vol_cor,2);

%% Visualise the new images %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
figure;
% Deformation Axial and Sagittal
subplot(131);
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
axis equal
axis off; grid off

% Deformation Axial and Coronal
subplot(132);

aa(:,:,1) = convert2u8(new_axial(:,:,slices(1)));
aa(:,:,2) = convert2u8(new_axial(:,:,slices(1)));
aa(:,:,3) = convert2u8(new_axial(:,:,slices(1)));

[X, Y, Z, triTexture] = compute_triTexture(axial_m{slices(1)},aa, [0 rows_ax-1],[0 cols_ax-1]);
surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal

cc(:,:,2) = convert2u8(new_coronal(:,:,slices(3)));
cc(:,:,1) = convert2u8(new_coronal(:,:,slices(3)));
cc(:,:,3) = convert2u8(new_coronal(:,:,slices(3)));

[X, Y, Z, triTexture] = compute_triTexture(cor_m{slices(3)},cc, [0 rows_cor-1],[0 cols_cor-1]);
surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal
axis off; grid off

% Deformation Sagittal and Coronal
subplot(133);

cc(:,:,2) = convert2u8(new_coronal(:,:,slices(3)));
cc(:,:,1) = convert2u8(new_coronal(:,:,slices(3)));
cc(:,:,3) = convert2u8(new_coronal(:,:,slices(3)));

[X, Y, Z, triTexture] = compute_triTexture(cor_m{slices(3)},cc, [0 rows_cor-1],[0 cols_cor-1]);
surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal


ss(:,:,2) = convert2u8(new_sagittal(:,:,slices(2)));
ss(:,:,1) = convert2u8(new_sagittal(:,:,slices(2)));
ss(:,:,3) = convert2u8(new_sagittal(:,:,slices(2)));

[X, Y, Z, triTexture] = compute_triTexture(sag_m{slices(2)},ss, [0 rows_sag-1],[0 cols_sag-1]);
surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal
axis off; grid off

%% Visualize Original %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Original Axial and Sagittal
figure;
subplot(131);
aa(:,:,1) = convert2u8(vol_ax(:,:,slices(1)));
aa(:,:,2) = convert2u8(vol_ax(:,:,slices(1)));
aa(:,:,3) = convert2u8(vol_ax(:,:,slices(1)));

[X, Y, Z, triTexture] = compute_triTexture(axial_m{slices(1)},aa,[0 size(views.axial,1)-1],[0 size(views.axial,2)-1]);
surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal

ss(:,:,2) = convert2u8(vol_sag(:,:,slices(2)));
ss(:,:,1) = convert2u8(vol_sag(:,:,slices(2)));
ss(:,:,3) = convert2u8(vol_sag(:,:,slices(2)));

[X, Y, Z, triTexture] = compute_triTexture(sag_m{slices(2)},ss,[0 size(views.sagittal,1)-1],[0 size(views.sagittal,2)-1]);
surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal
axis off; grid off


% Original Coronal and Sagittal
subplot(132);
cc(:,:,2) = convert2u8(vol_cor(:,:,slices(3)));
cc(:,:,1) = convert2u8(vol_cor(:,:,slices(3)));
cc(:,:,3) = convert2u8(vol_cor(:,:,slices(3)));

[X, Y, Z, triTexture] = compute_triTexture(cor_m{slices(3)},cc,[0 size(views.coronal,1)-1],[0 size(views.coronal,2)-1]);
surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal


ss(:,:,2) = convert2u8(vol_sag(:,:,slices(2)));
ss(:,:,1) = convert2u8(vol_sag(:,:,slices(2)));
ss(:,:,3) = convert2u8(vol_sag(:,:,slices(2)));

[X, Y, Z, triTexture] = compute_triTexture(sag_m{slices(2)},ss,[0 size(views.sagittal,1)-1],[0 size(views.sagittal,2)-1]);
surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal
axis off; grid off

% Original Coronal and Axial
subplot(133);
aa(:,:,1) = convert2u8(vol_ax(:,:,slices(1)));
aa(:,:,2) = convert2u8(vol_ax(:,:,slices(1)));
aa(:,:,3) = convert2u8(vol_ax(:,:,slices(1)));

[X, Y, Z, triTexture] = compute_triTexture(axial_m{slices(1)},aa,[0 size(views.axial,1)-1],[0 size(views.axial,2)-1]);
surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal

cc(:,:,2) = convert2u8(vol_cor(:,:,slices(3)));
cc(:,:,1) = convert2u8(vol_cor(:,:,slices(3)));
cc(:,:,3) = convert2u8(vol_cor(:,:,slices(3)));

[X, Y, Z, triTexture] = compute_triTexture(cor_m{slices(3)},cc,[0 size(views.coronal,1)-1],[0 size(views.coronal,2)-1]);
surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal
axis off; grid off
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Visualize the new images using the 3D viewer
images.axial      = new_axial;
images.axial_info = axial_m;

images.sagittal      = new_sagittal;
images.sagittal_info = sag_m;

images.coronal      = new_coronal;
images.coronal_info = cor_m;

view3d(images);


