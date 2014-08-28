%% Compute different gradients
global gradx_ax
global grady_ax
global gradx_sag
global grady_sag
global gradx_cor
global grady_cor

global vol_ax
global vol_sag
global vol_cor

total_ax = size(views.axial,3);
total_sag = size(views.sagittal,3);
total_cor = size(views.coronal,3);

% for i=1:total_ax
% 
%     [grady_ax(:,:,i) grady_ax(:,:,i)] = gradient(vol_ax(:,:,i));
% 
% end
% 
% for i=1:total_sag
%    
%     [grady_sag(:,:,i) gradx_sag(:,:,i)] = gradient(vol_sag(:,:,i));
% 
% end
% 
% for i=1:total_cor
%     
%     [grady_cor(:,:,i) gradx_cor(:,:,i)] = gradient(vol_cor(:,:,i));
%     
% end

sigm = 5;
G = fspecial('gaussian',[2*sigm+1 2*sigm+1],sigm);
[dx, dy] = gradient(fspecial('gauss',[2*sigm+1 2*sigm+1],sigm)); % G is a 2D gaussain

for i=1:total_ax

    grady_ax(:,:,i) = imfilter(vol_ax(:,:,i), dy,'same');
    gradx_ax(:,:,i) = imfilter(vol_ax(:,:,i), dx,'same');

end

for i=1:total_sag
   
    grady_sag(:,:,i) = imfilter(vol_sag(:,:,i), dy, 'same');
    gradx_sag(:,:,i) = imfilter(vol_sag(:,:,i), dx, 'same');

end
for i=1:total_cor
    
    grady_cor(:,:,i) = imfilter(vol_cor(:,:,i), dy, 'same');
    gradx_cor(:,:,i) = imfilter(vol_cor(:,:,i), dx, 'same');
end