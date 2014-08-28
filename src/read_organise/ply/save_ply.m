function save_ply(vertices, faces, combine_ply)

fileId = fopen(combine_ply, 'wt');

% writing the header
fprintf(fileId, 'ply\n');
fprintf(fileId, 'format ascii 1.0\n');
fprintf(fileId, 'comment created by safir\n');

% writing the elements
fprintf(fileId, 'element vertex %d\n', size(vertices,1));%var

fprintf(fileId, 'property float32 x\n');
fprintf(fileId, 'property float32 y\n');
fprintf(fileId, 'property float32 z\n');

fprintf(fileId, 'element face %d\n', size(faces,1));%t
fprintf(fileId, 'property list int32 int32 vertex_indices\n');

fprintf(fileId, 'property uchar red\n');
fprintf(fileId, 'property uchar green\n');
fprintf(fileId, 'property uchar blue\n');

fprintf(fileId, 'end_header\n');

% writing the vertexes

for i=1:size(vertices,1)
    fprintf(fileId, '%f %f %f \n', vertices(i,1),vertices(i,2),vertices(i,3));
end

%     % writing the faces
for i = 1 : size(faces,1)
    fprintf(fileId, '3 %d %d %d %u %u %u\n', int32(faces(i,1)), int32(faces(i,2)), int32(faces(i,3)), 0, 255, 255);
end
fclose(fileId);