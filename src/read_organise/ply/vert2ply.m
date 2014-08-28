function vert2ply(mat_file,ply_file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Given a .mat file containing 'vertex' as cells, write them in a txt file
%%  given by txt_file
%%
%%  NOTE: In the mat file the points are in the frame coordinates, and the
%%        ply file will have the points in the same coordinates system !
%%
%%  Inputs: 1. mat_file -> *.mat with 'vertex' as cells
%%          2. ply_file -> name of the ply file where 
%%                         to write the vertex info 
%% EXECUTE
%% vert2ply('resources/manual_seg/LASSALAS/lassalas_uterus_ax.mat','lassalas_uterus_ax.ply');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
a = load(mat_file);

vert = a.vertex;

%% number of slices
N = length(vert);

for i = 1:N
    N_slice(i) = size(vert{i},1); % # vertex per slice
end
N_total = sum(N_slice); % # total vertex

fileId = fopen(ply_file, 'wt');

% writing the header
fprintf(fileId, 'ply\n');
fprintf(fileId, 'format ascii 1.0\n');
fprintf(fileId, 'comment created by safir\n');

% writing the elements
fprintf(fileId, 'element vertex %d\n', N_total);
% [x y z] x,y image coordinates; z the current slice
fprintf(fileId, 'property float32 x\n');
fprintf(fileId, 'property float32 y\n');
fprintf(fileId, 'property float32 z\n');

fprintf(fileId, 'end_header\n');

% writing the vertexes
for i=1:N  
  for j=1:size(vert{i},1)
    fprintf(fileId, '%f %f %f \n', vert{i}(j,1), vert{i}(j,2),i);
  end
  
end

fclose(fileId);