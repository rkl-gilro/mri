function ASC = get_orientation(orientation)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Compute the image orientation: axial, sagittal, coronal, given the dicom
%% tag 'ImageOrientationPatient'
%%
%% Inputs:  - orientation -> 1X6 vector, it contains two vectors of 3 
%%                           components, direction cosine of the image frame with
%%                           respect the RCS.
%% Outputs: - ASC -> string with 3 possible different values
%%                   'a' axial, 's' sagittal, 'c' coronal
%%
%% Execute:
%% - [im, file] = dicom('resources/patients/patient1/exp0000', 0);
%% - asc = get_orientation( file.ImageOrientationPatient );
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

thr = 0.8; 

x = orientation(1:3);
y = orientation(4:end);

%% Get the major axis in the both directions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% x direction
if(x(1) < 0)
    x_orientation(1) = 'R';
else
    x_orientation(1) = 'L';
end
if(x(2) < 0)
    x_orientation(2) = 'A';
else
    x_orientation(2) = 'P';
end
if(x(3) < 0)
    x_orientation(3) = 'F';
else
    x_orientation(3) = 'H';
end

if(abs(x(1)) > thr)
    axis_x = x_orientation(1);
elseif(abs(x(2)) > thr)
    axis_x = x_orientation(2);
elseif(abs(x(3)) > thr)
    axis_x = x_orientation(3);
end
    
% y direction
if(y(1) < 0)
    y_orientation(1) = 'R';
else
    y_orientation(1) = 'L';
end
if(y(2) < 0)
    y_orientation(2) = 'A';
else
    y_orientation(2) = 'P';
end
if(y(3) < 0)
    y_orientation(3) = 'F';
else
    y_orientation(3) = 'H';
end

if(abs(y(1)) > thr)
    axis_y = y_orientation(1);
elseif(abs(y(2)) > thr)
    axis_y = y_orientation(2);
elseif(abs(y(3)) > thr)
    axis_y = y_orientation(3);
end

%% Compute the orientation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if( (strcmp(axis_x,'R') ||  strcmp(axis_x,'L'))  &&  (strcmp(axis_y,'P') || strcmp(axis_y,'A')) )
    ASC = 'A';
elseif( (strcmp(axis_y,'R') || strcmp(axis_y,'L')) && (strcmp(axis_x,'P') || strcmp(axis_x,'A')) )
    ASC = 'A';
    
elseif( (strcmp(axis_x,'R') || strcmp(axis_x,'L')) && (strcmp(axis_y,'F') || strcmp(axis_y,'H')) )
    ASC = 'C';
elseif( (strcmp(axis_y,'R') || strcmp(axis_y,'L')) && (strcmp(axis_x,'F') || strcmp(axis_x,'H')) )
    ASC = 'C';
    
elseif( (strcmp(axis_x,'A') || strcmp(axis_x,'P')) && (strcmp(axis_y,'F') || strcmp(axis_y,'H')) )
    ASC = 'S';
elseif( (strcmp(axis_y,'A') || strcmp(axis_y,'P')) && (strcmp(axis_x,'F') ||  strcmp(axis_x,'H')) )
    ASC = 'S';
end