function [X, Y, Z, triTexture] = compute_triTexture( patient, img, i, j)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Returns X, Y, Z and the texture for plotting the image in 
%%  the 3D space
%%
%%  Inputs:  1. patient  -> 4x3 matrix, the transformation from 2D to 3D from 
%%                         DICOM information
%%                       -> struct that contains the DICOM info to compute the    
%%                          4x3 matrix
%%           2.img -> matrix that contains the image, it is used for the texture
%%           3. i  -> 1x2 array of the rows image range, e.g [0 511]
%%           4. j  -> 1x2 array of the cols image range, e.g [0 511]
%%
%%  Outputs: 1. X  -> 4x1 array of the x components of the points
%%                    that define the image plane
%%           2. Y  -> idem 
%%           3. Z  -> idem 
%%           4. triTexture  -> the image transform in the 3D space
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

x = [];
y = [];
z = [];

%% Check if patient is a struct or a matrix %%
if isstruct(patient)
    
    M = zeros(4,3);
    
    M(1:3,3) = patient.ImagePositionPatient;
    M(4,3) = 1;
    M(1:3,1) = patient.ImageOrientationPatient(1:3).*patient.PixelSpacing(1);
    M(1:3,2) = patient.ImageOrientationPatient(4:6).*patient.PixelSpacing(2);
else
    M = patient;
end

%% Get the 3D locations of the corners of the image %%
for k=1:length(i)
    for l=1:length(j)
    p = M * [j(k)-1 i(l)-1 1]';% 0

        x = [x p(1)];
        y = [y p(2)];
        z = [z p(3)];

    end
end

[X, Y, Z, triTexture] = compute_RCS_img(x, y, z, img);

  
