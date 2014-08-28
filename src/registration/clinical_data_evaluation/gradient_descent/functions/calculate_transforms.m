function [M M1 X Y Z plane] = calculate_transforms(view_info, view_size, range1, range2, ortho, view)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Calculates the transformation from 2D image coordinates to 3D space RCS
%%
%%  Inputs:  1. view_info -> cell containing the DICOM info for each slice
%%                           in one view
%%           2. view_size -> [m n s] is the size of the volume
%%           3. range1    -> [c1 c2] min and max value of columns
%%           4. range2    -> [r1 r2] min and max value of rows
%%           5. ortho     -> 0/1 -- not orthogonal/orthogonal
%%           6. view      -> scalar 1-axial 2-sagittal 3-coronal
%%
%%  Outputs: 1. M     -> cell that contains 's' 4x3 matrix from 2D (i,j) to 3D (x,y,z)
%%           2. M1    -> cell that contains 's' 3x4 matrix inverse of M
%%           3. X     -> 4xM matrix contains the x coordinates of the four 
%%                      points that defines the image plane in 3D
%%           4. Y     -> idem, y coordinates
%%           5. Z     -> idem, y coordinates
%%           6. plane -> sx4 matrix which each row contains
%%                      (A,B,C,D) of the plane Ax + By + Cz + D =0
%%
%% NOTE: in case we crop the volume, range1 & range2 will be different than
%%       [0 cols-1] & [0 rows-1]
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Initialize variables
X = [];
Y = [];
Z = [];

M  = cell(1,view_size(3));
M1 = cell(1,view_size(3));

plane = zeros(view_size(3),4);

for i = 1:view_size(3)
    
    [M{i}, M1{i}] = compute_M_M1(view_info{i}, ortho, view);
    
    [x,y,z] = calculate4corners( M{i}, range1, range2 ); % [0 view_size(2)-1], [0 view_size(1)-1]
    
    X = [X x'];
    Y = [Y y'];
    Z = [Z z'];
    
    if i == 1
        
        N = cross( -[X(1,i) Y(1,i) Z(1,i)] + [X(2,i) Y(2,i) Z(2,i)], -[X(1,i) Y(1,i) Z(1,i)] + [X(3,i) Y(3,i) Z(3,i)]);
        N = N./norm(N);
        
    end
    
    plane(i,:) = [N -(N(1)*X(4,i) + N(2)*Y(4,i) + N(3)*Z(4,i))]; % (A,B,C,D) of the plane Ax + By + Cz + D =0
    
end