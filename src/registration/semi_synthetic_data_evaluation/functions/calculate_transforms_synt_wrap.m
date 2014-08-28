function [M M1 X Y Z plane] = calculate_transforms_synt_wrap(view_info, view_size, range1, range2, range3, ortho, M_ax, view)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Calculates the transformation from 2D image coordinates to 3D space RCS
%%  Using synthetic data
%%
%%  Inputs:  1. view_info -> struct contains the DICOM info for each slice
%%                           in one view
%%           2. view_size -> [m n s] is the size of the volume
%%           3. range1    -> [c1 c2] min and max value of columns
%%           4. range2    -> [r1 r2] min and max value of rows
%%           5. range3    -> [s1 s2] min and max value of slices
%%           6. ortho     -> 0/1 -- not orthogonal/orthogonal
%%           7. M_ax      -> transformation matrix from the axial view 
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

iter = 0;

for i = range3(1):range3(2)
    
    if view == 2
        tmp = M_ax * [i-1 0 1]';
    else
        tmp = M_ax * [0 i-1 1]';
    end
    
    iter = iter + 1;
    [M{iter}, M1{iter}, ~] = compute_M_M1_synt(view_info, ortho, view, tmp(1:3));
    
    [x,y,z] = calculate4corners( M{iter}, range1, range2 ); %[0 view_size(2)-1], [0 view_size(1)-1]
    
    X = [X x'];
    Y = [Y y'];
    Z = [Z z'];
    
    if iter == 1
        
        N = cross( -[X(1,iter) Y(1,iter) Z(1,iter)] + [X(2,iter) Y(2,iter) Z(2,iter)], -[X(1,iter) Y(1,iter) Z(1,iter)] + [X(3,iter) Y(3,iter) Z(3,iter)]);
        N = N./norm(N);
        
    end
    
    plane(iter,:) = [N -(N(1)*X(4,iter) + N(2)*Y(4,iter) + N(3)*Z(4,iter))]; % (A,B,C,D) of the plane Ax + By + Cz + D =0
    
end