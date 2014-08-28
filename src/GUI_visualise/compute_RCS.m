function [X, Y, Z, triTexture] = compute_RCS(patient,img)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = zeros(4,4);

M(1:3,4) = patient.ImagePositionPatient;
M(4,4) = 1;
M(1:3,1) = patient.ImageOrientationPatient(1:3).*patient.PixelSpacing(1);
M(1:3,2) = patient.ImageOrientationPatient(4:6).*patient.PixelSpacing(2);

i = [0 511];
j = [0 511];

x = [];
y = [];
z = [];


for k=1:length(i)
    for l=1:length(j)
    p = M*[j(k) i(l) 0 1]';

        x = [x p(1)];
        y = [y p(2)];
        z = [z p(3)];

    end
end

P1 = [x(1) y(1) z(1)];
P2 = [x(2) y(2) z(2)];
P3 = [x(3) y(3) z(3)];
P4 = [x(4) y(4) z(4)];

x = [P1(1) P4(1) P2(1)];   %# [xorigin xA xB] coordinates in 3-D space
y = [P1(2) P4(3) P2(2)];   %# [yorigin yA yB] coordinates in 3-D space
z = [P1(3) P4(3) P2(3)];   %# [zorigin zA zB] coordinates in 3-D space

origin = [512 1];  %# Vertex of triangle in image space
U = [0 511];       %# Vector from origin to point A in image space
V = [-511 511]; 

A = origin+U;  %# Point A
B = origin+V;  %# Point B
C = B-U;     %# Point C

[nRows,nCols,nPages] = size(img);  %# Image dimensions
inputCorners = [origin; ...        %# Corner coordinates of input space
                A; ...
                B; ...
                C];
            
outputCorners = [origin; ...      %# Corner coordinates of output space
                 A; ...
                 B; ...
                 C];            
             
tform = maketform('projective',...  %# Make the transformation structure
                  inputCorners,...
                  outputCorners);
triTexture = imtransform(img,tform,'bicubic',...  %# Transform the image
                         'xdata',[1 nCols],...
                         'ydata',[1 nRows],...
                         'size',[nRows nCols]);
x = [P3(1) P4(1) P2(1) P1(1)];   %# [xorigin xA xB] coordinates in 3-D space
y = [P3(2) P4(2) P2(2) P1(2)];   %# [yorigin yA yB] coordinates in 3-D space
z = [P3(3) P4(3) P2(3) P1(3)];   %# [zorigin zA zB] coordinates in 3-D space

index = [4 1;3 2]; %# Index used to create 2-by-2 surface coordinates
X = x(index);      %# x coordinates of surface
Y = y(index);       %# y coordinates of surface
Z = z(index);        %# z coordinates of surface