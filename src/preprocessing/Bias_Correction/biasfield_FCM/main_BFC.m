% Compile the c-code
% mex src/preprocessing/Bias_Correction/biasfield_FCM/BCFCM3D.c -v;
% Load test image
Y = single(views.sagittal(:,:,:));

% Class prototypes (means)
v = [10;100;500;1000];
% Do the fuzzy clustering
[B,U]=BCFCM3D(Y,v,struct('maxit',25,'epsilon',1e-8));
% Show results
% figure,
% subplot(2,2,1), imshow(Y,[]), title('input image');
% subplot(2,2,2), imshow(U,[]), title('Partition matrix');
% subplot(2,2,3), imshow(B,[]), title('Estimated biasfield');
% subplot(2,2,4), imshow(Y-B,[]), title('Corrected image');
