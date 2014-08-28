%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Main Laws Filter
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Select input image
im = views.axial(:,:,10);
%im = views.sagittal(:,:,10);
%im = views.coronal(:,:,12);

%% Compute the masks
%masks = Laws_filter_definition;

%% Select a mask
mask = masks.S3S3;

%% Calculate the filtered image
statistics = Laws_filter(im,mask);

%% Show results
figure;
subplot(121);imshow(statistics.mean,[]);title('Mean');
subplot(122);imshow(statistics.stdev,[]);title('Standard Deviation');