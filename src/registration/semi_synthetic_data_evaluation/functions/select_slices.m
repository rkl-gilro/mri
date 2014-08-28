%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Select slices compute the slice selection once the deformation has been
%% applied to each of the volumes
%%
%%  1. After it is also computed the gradient of the new slices.
%%  2. The default number of slices per view is 25.
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global X_ax_v
global Y_ax_v
global Z_ax_v

global X_sag_v
global Y_sag_v
global Z_sag_v

global X_cor_v
global Y_cor_v
global Z_cor_v

global X_ax_def
global Y_ax_def
global Z_ax_def

global X_sag_def
global Y_sag_def
global Z_sag_def

global X_cor_def
global Y_cor_def
global Z_cor_def

global lava_axM
global lava_axM_1
global lava_sagM
global lava_sagM_1
global lava_corM
global lava_corM_1

global axial_M
global axial_M1
global sag_M
global sag_M1
global cor_M
global cor_M1

global slices_ax
global slices_sag
global slices_cor

global vol_ax_eval
global vol_sag_eval
global vol_cor_eval

global vol_ax_eval1
global vol_sag_eval1
global vol_cor_eval1

%% Select a few slices, like in T2W scans

global size1 size2 size3

n_slices_ax = 25;
slices_ax = 5:round(size3/n_slices_ax):size3-4;

st_sag   = 100;
end_sag = size2 - st_sag;
n_slices_sag = 25;
slices_sag = st_sag:ceil((end_sag-st_sag)/n_slices_sag):end_sag;

st_cor   = 60;
end_cor = size1 - 100;%st_cor;
n_slices_cor = 25;
slices_cor = st_cor:ceil((end_cor-st_cor)/n_slices_cor):end_cor;

ax  = slices_ax;
sag = slices_sag;
cor = slices_cor;

disp('--------- Select slices for axial -----')
for i = 1:length(slices_ax)
    ind = slices_ax(i);

    X_ax_def = [X_ax_def X_ax_v(:,ind)];
    Y_ax_def = [Y_ax_def Y_ax_v(:,ind)];
    Z_ax_def = [Z_ax_def Z_ax_v(:,ind)];
    
    vol_ax_eval(:,:,i) = vol_ax_eval1(:,:,ind); 
    
    axial_M{i}  = lava_axM{ind}; 
    axial_M1{i} = lava_axM_1{ind}; 

end


disp('--------- Select slices for sagittal -----')
for i = 1:length(slices_sag)
    ind = slices_sag(i);
    
    X_sag_def = [X_sag_def X_sag_v(:,ind)];
    Y_sag_def = [Y_sag_def Y_sag_v(:,ind)];
    Z_sag_def = [Z_sag_def Z_sag_v(:,ind)];
    
    vol_sag_eval(:,:,i) = vol_sag_eval1(:,:,ind); 
    
    sag_M{i}  = lava_sagM{ind}; 
    sag_M1{i} = lava_sagM_1{ind}; 

end


disp('--------- Select slices for coronal -----')
for i = 1:length(slices_cor)
    ind = slices_cor(i);

    X_cor_def = [X_cor_def X_cor_v(:,ind)];
    Y_cor_def = [Y_cor_def Y_cor_v(:,ind)];
    Z_cor_def = [Z_cor_def Z_cor_v(:,ind)];
    
    vol_cor_eval(:,:,i) = vol_cor_eval1(:,:,ind); 
    
    cor_M{i}  =  lava_corM{ind}; 
    cor_M1{i} = lava_corM_1{ind}; 

end

global grady_ax
global gradx_ax

global grady_sag
global gradx_sag

global grady_cor
global gradx_cor

disp('--------- Calculate image gradients -----')
%% Image gradients %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
H = fspecial('prewitt');

sigm = 7;
G = fspecial('gaussian',[2*sigm+1 2*sigm+1],sigm);
[dx, dy] = gradient(fspecial('gauss',[2*sigm+1 2*sigm+1],sigm)); % G is a 2D gaussain

        
for i=1:size(vol_ax_eval,3)
    
    grady_ax(:,:,i) = imfilter((vol_ax_eval(:,:,i)), dy,'same'); % conv2(vol_ax(:,:,i), H', 'same');%
    gradx_ax(:,:,i) = imfilter((vol_ax_eval(:,:,i)), dx,'same');

end

for i=1:size(vol_sag_eval,3)
    
    grady_sag(:,:,i) = imfilter((vol_sag_eval(:,:,i)), dy, 'same');
    gradx_sag(:,:,i) = imfilter((vol_sag_eval(:,:,i)), dx, 'same');

end
for i=1:size(vol_cor_eval,3)
    
    grady_cor(:,:,i) = imfilter((vol_cor_eval(:,:,i)), dy, 'same');
    gradx_cor(:,:,i) = imfilter((vol_cor_eval(:,:,i)), dx, 'same');
    
end