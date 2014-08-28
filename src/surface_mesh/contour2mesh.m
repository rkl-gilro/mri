function [faces, vertices] = contour2mesh(input_vert, patient, view, ply_file, show)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  From a given set of contours (all the frames from one view) in a mat file
%%  from a given patient, compute the surface triangulation and save it in
%%  a ply file
%%
%%  NOTE: It is important to note that in the mat_file the contours are in 
%%        frame coordinates, rather than in the output ply file, the contour
%%        points are given in the RCS !
%%
%%  Inputs: 1. input_vert -> 1. (string) of the path of the .mat file that contains
%%                         the structure 'vertex', a cell array where each
%%                         cell contains the vertex of the contour in one slice
%%                    or     2. (cell) a cell array where each
%%                         cell contains the vertex of the contour in one slice
%%          2. patient ->  (cell) with the DICOM info for each slice
%%          3. view -> 1-2 3 := axial-sagittal-coronal
%%          4. ply_file(optional) -> (string) path of the ply file to save the 
%%                                   vertex and faces information
%%
%%  Output: 1. faces -> Nfx3 matrix, where each row is a face and the columns
%%                              correspond to the vertex of this face
%%          2. vertices -> Nv x 3 matrix with rows number of contour points, 
%%                    and columns (x,y,z) coordinates
%%
%% EXECUTE
%% main_read;
%% [faces, vertices] = contour2mesh('resources/manual_seg/LASSALAS/lassalas_uterurs_sag2.mat',...
%%                         views.sagittal_info,'lassalas_sag.ply');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

save_file = 1;

%% Check the number of inputs
if nargin < 4
    show      = 0;
    save_file = 0;
end

if nargin < 3
    save_file = 0;
    show      = 0;
end

%% If input_vert is a cell, the vertex are directly the ones provided by
%% input_vert
if iscell(input_vert)
    vert = input_vert;
end
%% If input_vert is a string, the file should be loaded to get the verices
if ischar(input_vert)
    a = load(input_vert);
    vert = a.vertex;
end

%% Initialize
x = [];
y = [];
z = [];

true_x = [];
true_y = [];
true_z = [];

%% Transform the contour points into the reference coordinate system (RCS)
for i = 1:length(patient)
    
    if ~isempty(vert{i})
        
        M = zeros(4,4);
        
        %% Compute transformation from pixel coordinates to RCS
        M(1:3,4) = patient{i}.ImagePositionPatient;
        M(4,4) = 1;
        M(1:3,1) = patient{i}.ImageOrientationPatient(1:3).*patient{i}.PixelSpacing(1);
        M(1:3,2) = patient{i}.ImageOrientationPatient(4:6).*patient{i}.PixelSpacing(2);
        
        %% Get the points in RCS
        p = M * [vert{i}(:,2) vert{i}(:,1) zeros(size(vert{i},1),1) ones(size(vert{i},1),1)]';
        
        true_x = [true_x p(1,:)];
        true_y = [true_y p(2,:)];
        true_z = [true_z p(3,:)];
        
%% Uncomment if you want to interpolate to get more points
        [Vertices, ~] = contour(vert{i});
        for k=1:size(Vertices,1)
            
            p = M * [Vertices(k,2) Vertices(k,1) 0 1]';
            
            x = [x p(1)];
            y = [y p(2)];
            z = [z p(3)];
            
        end
    end
    
end

var(:,1) = x(:)'; % true_x(:)'
var(:,2) = y(:)'; % true_y(:)'
var(:,3) = z(:)'; % true_z(:)'

%% Use 'unique' to discard the repeted points
[~,I,~] = unique(var,'first','rows');
I = sort(I);
var = var(I,:);

%% Calculate the surface based on 'PowerCrust Algorithm'
if view == 2
    [t] = Crust(var);
else 
    [t] = MyCrust(var);
end

faces = t;
vertices = var;


%% Plot of the output triangulation
if show
    figure;
    hold on
    title('Output Triangulation','fontsize',14)
    axis equal
    trisurf(t,var(:,1),var(:,2),var(:,3),'facecolor','c','edgecolor','b')%plot della superficie trattata
    hold on
end

%% Save the contours in a ply file 
if save_file
    
    save_ply(var, t, ply_file);
    
end
