% % Compile the c-code
  mex src/3Dsegmentation/FCM/BCFCM3D.c -v
  
% Load test volume
  D = read_dbt('resources/651_slices1mm_order/',.2,0);
  
% Convert to single
  D = im2single(D);

% Class prototypes (means)
  v = [0 0.12 0.18];
  
% Do the fuzzy clustering
  [B,U] = BCFCM3D(D,v,struct('maxit',5,'epsilon',1e-5,'sigma',1));
  
  cl_vol = classes2gray(U);

% Show results
  figure, 
  subplot(2,2,1), imshow(D(:,:,15),[]), title('A slice of input volume');
  subplot(2,2,2), imshow(squeeze(U(:,:,15,1:3))), title('3 classes of the partition matrix');
  subplot(2,2,3), imshow(B(:,:,15),[]), title('Estimated biasfield');
  subplot(2,2,4), imshow(D(:,:,15)-B(:,:,15),[]), title('Corrected slice');