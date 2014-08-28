function contour2ply(mat_file, patient, ply_file)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

a = load(mat_file);

vert = a.vertex;

%% number of slices
N = length(vert);

x = [];
y = [];
z = [];

for i = 1:length(patient)
    
    if ~isempty(vert{i})
        M = zeros(4,4);
        
        M(1:3,4) = patient{i}.ImagePositionPatient;
        M(4,4) = 1;
        M(1:3,1) = patient{i}.ImageOrientationPatient(1:3).*patient{i}.PixelSpacing(1);
        M(1:3,2) = patient{i}.ImageOrientationPatient(4:6).*patient{i}.PixelSpacing(2);
        

        for k=1:size(vert{i},1)
            
            p = M * [vert{i}(k,1) vert{i}(k,2) 0 1]';
            
            x = [x p(1)];
            y = [y p(2)];
            z = [z p(3)];
            
        end
    end
    
end

fileId = fopen(ply_file, 'wt');

% writing the header
fprintf(fileId, 'ply\n');
fprintf(fileId, 'format ascii 1.0\n');
fprintf(fileId, 'comment created by safir\n');

% writing the elements
fprintf(fileId, 'element vertex %d\n', length(x));
% [x y z] x,y image coordinates; z the current slice
fprintf(fileId, 'property float32 x\n');
fprintf(fileId, 'property float32 y\n');
fprintf(fileId, 'property float32 z\n');

fprintf(fileId, 'end_header\n');

% writing the vertexes

for i=1:length(x)
    
    fprintf(fileId, '%f %f %f \n', x(i), y(i),z(i));
    
end

fclose(fileId);