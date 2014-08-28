%%  maxdiff_intensity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Compute the maximum difference between intersection before the 
%%  registration is done
%%
%%  NOTE: load the registration
%% 1. load('/home/raquel/Documents/repositories/mri-segmentation/resources/registration/lassalas/5gradt10_74_.01/lassalas.mat')
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


measure_errors.total_diff_a;
measure_errors.total_diff_b;

[~,~,t] = size(measure_errors.diff_ax1_b); 

int1 = sum(measure_errors.diff_ax1_b,3);
int2 = sum(measure_errors.diff_cor2_b,3);
int3 = sum(measure_errors.diff_sag3_b,3);

max_val_b = max([max(vol_ax(:)) max(vol_sag(:)) max(vol_cor(:))]);
%% Compute the maximum intersection slice between axial&sagittal
[v1 ind1] = max(int1(:));
% [rind1_1 rind2_1] = ind2sub(size(int1),ind1);
rind1_1 = 5;
rind2_1 = 18;
disp(['Before: axial & sagittal ', ...
num2str(mean(squeeze(abs(measure_errors.diff_ax1_b(rind1_1,rind2_1,:)-measure_errors.diff_sag1_b(rind1_1,rind2_1,:))))), '+', ...
num2str(std(squeeze(abs(measure_errors.diff_ax1_b(rind1_1,rind2_1,:)-measure_errors.diff_sag1_b(rind1_1,rind2_1,:))))) ])

%% Compute the maximum intersection slice between axial&coronal
[v2 ind2] = max(int2(:));
[rind1_2 rind2_2] = ind2sub(size(int2),ind2);

disp(['Before: axial & coronal ', ...
num2str(mean(squeeze(abs(measure_errors.diff_ax2_b(rind1_2,rind2_2,:) - measure_errors.diff_cor2_b(rind1_2,rind2_2,:))))), '+', ...
num2str(std(squeeze(abs(measure_errors.diff_ax2_b(rind1_2,rind2_2,:) - measure_errors.diff_cor2_b(rind1_2,rind2_2,:))))) ])

%% Compute the maximum intersection slice between coronal&sagittal
[v3 ind3] = max(int3(:));
[rind1_3 rind2_3] = ind2sub(size(int3),ind3);

disp(['Before: coronal & sagittal ', ...
num2str(mean(squeeze(abs(measure_errors.diff_sag3_b(rind1_3,rind2_3,:) - measure_errors.diff_cor3_b(rind1_3,rind2_3,:))))), '+', ...
num2str(std(squeeze(abs(measure_errors.diff_sag3_b(rind1_3,rind2_3,:) - measure_errors.diff_cor3_b(rind1_3,rind2_3,:))))) ])

figure;plot(0:t-1, squeeze(measure_errors.diff_ax1_b(rind1_1,rind2_1,:)),'b');hold on
plot(0:t-1, squeeze(measure_errors.diff_sag1_b(rind1_1,rind2_1,:)),'r');hold on
axis([0, t-1, 0, max_val_b])
legend('Axial','Sagittal')
xlabel('Intersections') % x-axis label
ylabel('Intensity values') % y-axis label

figure;plot(0:t-1, squeeze(measure_errors.diff_ax2_b(rind1_2,rind2_2,:)),'b');hold on
plot(0:t-1, squeeze(measure_errors.diff_cor2_b(rind1_2,rind2_2,:)),'r');hold on
axis([0, t-1, 0, max_val_b])
legend('Axial','Coronal')
xlabel('Intersections') % x-axis label
ylabel('Intensity values') % y-axis label

figure;plot(0:t-1, squeeze(measure_errors.diff_sag3_b(rind2_3,rind2_3,:)),'b');hold on
plot(0:t-1, squeeze(measure_errors.diff_cor3_b(rind2_3,rind1_3,:)),'r');hold on
axis([0, t-1, 0, max_val_b])
legend('Sagittal','Coronal')
xlabel('Intersections') % x-axis label
ylabel('Intensity values') % y-axis label

[rows_ax, cols_ax,~]   = size(vol_ax);
[rows_sag, cols_sag,~] = size(vol_sag);
[rows_cor, cols_cor,~] = size(vol_cor);

%% Axial & Sagittal
figure;
slices = [rind1_1, rind2_1, 1]
aa(:,:,1) = convert2u8(vol_ax(:,:,slices(1)));
aa(:,:,2) = convert2u8(vol_ax(:,:,slices(1)));
aa(:,:,3) = convert2u8(vol_ax(:,:,slices(1)));
[X, Y, Z, triTexture] = compute_triTexture(axial_m{slices(1)},aa, [0 rows_ax-1],[0 cols_ax-1]);
surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal

ss(:,:,2) = convert2u8(vol_sag(:,:,slices(2)));
ss(:,:,1) = convert2u8(vol_sag(:,:,slices(2)));
ss(:,:,3) = convert2u8(vol_sag(:,:,slices(2)));
[X, Y, Z, triTexture] = compute_triTexture(sag_m{slices(2)},ss, [0 rows_sag-1],[0 cols_sag-1]);
surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal
axis off

%% Axial & Coronal
figure;
slices = [rind1_2, 1, rind2_2]
aa(:,:,1) = convert2u8(vol_ax(:,:,slices(1)));
aa(:,:,2) = convert2u8(vol_ax(:,:,slices(1)));
aa(:,:,3) = convert2u8(vol_ax(:,:,slices(1)));
[X, Y, Z, triTexture] = compute_triTexture(axial_m{slices(1)},aa, [0 rows_ax-1],[0 cols_ax-1]);
surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal

ss(:,:,2) = convert2u8(vol_cor(:,:,slices(3)));
ss(:,:,1) = convert2u8(vol_cor(:,:,slices(3)));
ss(:,:,3) = convert2u8(vol_cor(:,:,slices(3)));
[X, Y, Z, triTexture] = compute_triTexture(cor_m{slices(3)},ss, [0 rows_cor-1],[0 cols_cor-1]);
surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal
axis off

%% Coronal & Sagittal
figure;
slices = [1, rind1_3, rind2_3]
aa(:,:,1) = convert2u8(vol_cor(:,:,slices(3)));
aa(:,:,2) = convert2u8(vol_cor(:,:,slices(3)));
aa(:,:,3) = convert2u8(vol_cor(:,:,slices(3)));
[X, Y, Z, triTexture] = compute_triTexture(cor_m{slices(3)},aa, [0 rows_cor-1],[0 cols_cor-1]);
surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal

ss(:,:,2) = convert2u8(vol_sag(:,:,slices(2)));
ss(:,:,1) = convert2u8(vol_sag(:,:,slices(2)));
ss(:,:,3) = convert2u8(vol_sag(:,:,slices(2)));
[X, Y, Z, triTexture] = compute_triTexture(sag_m{slices(2)},ss, [0 rows_sag-1],[0 cols_sag-1]);
surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal
axis off


%%  maxdiff_intensity
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Compute the maximum intensity difference intersection after registration
%%  Actually it is used the value computed before registration for comparing
%%  the results
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


measure_errors.total_diff_a;
measure_errors.total_diff_b;

[~,~,t] = size(measure_errors.diff_ax1_a); 

int1 = sum(measure_errors.diff_ax1_a,3);
int2 = sum(measure_errors.diff_cor2_a,3);
int3 = sum(measure_errors.diff_sag3_a,3);

%% Compute the maximum intersection slice between axial&sagittal
% [v1 ind1] = max(int1(:));
% [rind1_1 rind2_1] = ind2sub(size(int1),ind1);

max_val_a = max([max(new_axial(:)) max(new_sagittal(:)) max(new_coronal(:))]);

disp(['After: axial & sagittal ', ...
num2str(mean(squeeze(abs(measure_errors.diff_ax1_a(rind1_1,rind2_1,:)-measure_errors.diff_sag1_a(rind1_1,rind2_1,:))))),'+', ...
num2str(std(squeeze(abs(measure_errors.diff_ax1_a(rind1_1,rind2_1,:)-measure_errors.diff_sag1_a(rind1_1,rind2_1,:))))) ])

%% Compute the maximum intersection slice between axial&coronal
% [v2 ind2] = max(int2(:));
% [rind1_2 rind2_2] = ind2sub(size(int2),ind2);

disp(['After: axial & coronal ', ...
num2str(mean(squeeze(abs(measure_errors.diff_ax2_a(rind1_2,rind2_2,:)-measure_errors.diff_cor2_a(rind1_2,rind2_2,:))))),'+', ...
num2str(std(squeeze(abs(measure_errors.diff_ax2_a(rind1_2,rind2_2,:)-measure_errors.diff_cor2_a(rind1_2,rind2_2,:))))) ])

%% Compute the maximum intersection slice between coronal&sagittal
% [v3 ind3] = max(int3(:));
% [rind1_3 rind2_3] = ind2sub(size(int3),ind3);

disp(['After: coronal & sagittal ', ...
num2str(mean(squeeze(abs(measure_errors.diff_sag3_a(rind1_3,rind2_3,:)-measure_errors.diff_cor3_a(rind1_3,rind2_3,:))))) ,'+', ...
num2str(std(squeeze(abs(measure_errors.diff_sag3_a(rind1_3,rind2_3,:)-measure_errors.diff_cor3_a(rind1_3,rind2_3,:))))) ])


figure;plot(0:t-1, squeeze(measure_errors.diff_ax1_a(rind1_1,rind2_1,:)),'b');hold on
plot(0:t-1, squeeze(measure_errors.diff_sag1_a(rind1_1,rind2_1,:)),'r');hold on
axis([0, t-1, 0, max_val_b])
legend('Axial','Sagittal')
xlabel('Intersections') % x-axis label
ylabel('Intensity values') % y-axis label

figure;plot(0:t-1, squeeze(measure_errors.diff_ax2_a(rind1_2,rind2_2,:)),'b');hold on
plot(0:t-1, squeeze(measure_errors.diff_cor2_a(rind1_2,rind2_2,:)),'r');hold on
axis([0, t-1, 0, max_val_b])
legend('Axial','Coronal')
xlabel('Intersections') % x-axis label
ylabel('Intensity values') % y-axis label

figure;plot(0:t-1, squeeze(measure_errors.diff_sag3_a(rind2_3,rind2_3,:)),'b');hold on
plot(0:t-1, squeeze(measure_errors.diff_cor3_a(rind2_3,rind1_3,:)),'r');hold on
axis([0, t-1, 0, max_val_b])
legend('Coronal','Sagittal')
xlabel('Intersections') % x-axis label
ylabel('Intensity values') % y-axis label

[rows_ax, cols_ax,~]   = size(new_axial);
[rows_sag, cols_sag,~] = size(new_sagittal);
[rows_cor, cols_cor,~] = size(new_coronal);

%% Axial & Sagittal
figure;
slices = [rind1_1, rind2_1, 1];
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
axis off

%% Axial & Coronal
figure;
slices = [rind1_2, 1, rind2_2];
aa(:,:,1) = convert2u8(new_axial(:,:,slices(1)));
aa(:,:,2) = convert2u8(new_axial(:,:,slices(1)));
aa(:,:,3) = convert2u8(new_axial(:,:,slices(1)));
[X, Y, Z, triTexture] = compute_triTexture(axial_m{slices(1)},aa, [0 rows_ax-1],[0 cols_ax-1]);
surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal

ss(:,:,2) = convert2u8(new_coronal(:,:,slices(3)));
ss(:,:,1) = convert2u8(new_coronal(:,:,slices(3)));
ss(:,:,3) = convert2u8(new_coronal(:,:,slices(3)));
[X, Y, Z, triTexture] = compute_triTexture(cor_m{slices(3)},ss, [0 rows_cor-1],[0 cols_cor-1]);
surf(X,Y,Z,triTexture,'FaceColor','texturemap','EdgeColor','none');
hold on
axis equal
axis off

%% Coronal & Sagittal
figure;
slices = [1, rind1_3, rind2_3];
aa(:,:,1) = convert2u8(new_coronal(:,:,slices(3)));
aa(:,:,2) = convert2u8(new_coronal(:,:,slices(3)));
aa(:,:,3) = convert2u8(new_coronal(:,:,slices(3)));
[X, Y, Z, triTexture] = compute_triTexture(cor_m{slices(3)},aa, [0 rows_cor-1],[0 cols_cor-1]);
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
axis off

%%%%%%%%%%%%%%
%% Plot bars
b1 = squeeze(abs(measure_errors.diff_ax1_b(rind1_1,rind2_1,:)-measure_errors.diff_sag1_b(rind1_1,rind2_1,:)));
a1 = squeeze(abs(measure_errors.diff_ax1_a(rind1_1,rind2_1,:)-measure_errors.diff_sag1_a(rind1_1,rind2_1,:)));

b2 = squeeze(abs(measure_errors.diff_ax2_b(rind1_2,rind2_2,:)-measure_errors.diff_cor2_b(rind1_2,rind2_2,:)));
a2 = squeeze(abs(measure_errors.diff_ax2_a(rind1_2,rind2_2,:)-measure_errors.diff_cor2_a(rind1_2,rind2_2,:)));

b3 = squeeze(abs(measure_errors.diff_sag3_b(rind1_3,rind2_3,:)-measure_errors.diff_cor3_b(rind1_3,rind2_3,:)));
a3 = squeeze(abs(measure_errors.diff_sag3_a(rind1_3,rind2_3,:)-measure_errors.diff_cor3_a(rind1_3,rind2_3,:)));

C    = cell(t*6,1);

C(1:3*t) = {'Before'};
C(3*t+1:6*t) = {'After'};

figure;
% myboxplot([ones(1,t) 2.*ones(1,t) 3.*ones(1,t) 4.*ones(1,t) 5.*ones(1,t) 6.*ones(1,t)]', ...
%            [b1' a1' b2' a2' b3' a3']', 'split',C,'style_tukey','leglocation','north' );
myboxplot([ones(1,3*t) 2.*ones(1,3*t) ]', ...
           [b1' a1' b2' a2' b3' a3']', 'split',C,'style_tukey','leg','auto' );
% boxplot([b1 a1 b2 a2 b3 a3],{'Before' 'After' 'Before' 'After' 'Before' 'After'});
