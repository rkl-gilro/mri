function plot_views(views,ax,sag,cor)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Displays the 3 planes in the reference coordinate system 
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


aa(:,:,1) = convert2u8(views.axial(:,:,ax));
aa(:,:,2) = zeros(size(views.axial(:,:,ax)));
aa(:,:,3) = zeros(size(views.axial(:,:,ax)));


ss(:,:,2) = convert2u8(views.sagittal(:,:,sag));
ss(:,:,1) = zeros(size(views.sagittal(:,:,sag)));
ss(:,:,3) = zeros(size(views.sagittal(:,:,sag)));

cc(:,:,3) = convert2u8(views.coronal(:,:,cor));
cc(:,:,1) = zeros(size(views.coronal(:,:,cor)));
cc(:,:,2) = zeros(size(views.coronal(:,:,cor)));



%% Axial %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

M = zeros(4,4);

M(1:3,4) = views.axial_info{ax}.ImagePositionPatient;
M(4,4) = 1;
M(1:3,1) = views.axial_info{ax}.ImageOrientationPatient(1:3).*views.axial_info{ax}.PixelSpacing(1);
M(1:3,2) = views.axial_info{ax}.ImageOrientationPatient(4:6).*views.axial_info{ax}.PixelSpacing(2);

% dot(views.axial_info{ax}.ImageOrientationPatient(1:3),views.axial_info{ax}.ImageOrientationPatient(4:6))
% norm(views.axial_info{ax}.ImageOrientationPatient(4:6))
% norm(views.axial_info{ax}.ImageOrientationPatient(1:3))


i = [0 511];
j = [0 511];

x = [];
y = [];
z = [];

figure;
for k=1:length(i)
    for l=1:length(j)
    p = M*[i(k) j(l) 0 1]';

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
img = aa;  %# Sample image for texture map

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
% index = [3 4; 2 1];  %# Index used to create 2-by-2 surface coordinates
%index = [1 4;2 3];
index = [4 1;3 2];
X1 = x(index);      %# x coordinates of surface
Y1 = y(index);       %# y coordinates of surface
Z1 = z(index);        %# z coordinates of surface

plot3(P1(1),P1(2),P1(3),'r*','MarkerSize',12);hold on
plot3(P2(1),P2(2),P2(3),'b*','MarkerSize',12);hold on
plot3(P3(1),P3(2),P3(3),'g*','MarkerSize',12);hold on
plot3(P4(1),P4(2),P4(3),'k*','MarkerSize',12);hold on

hSurface = surf(X1,Y1,Z1,triTexture,...          %# Plot texture-mapped surface
                'FaceColor','texturemap',...
                'EdgeColor','none');
axis equal 

% normal = cross(P1-P2, P1-P3);
% D = -(normal(1)*P4(1) + normal(2)*P4(2)  + normal(3)*P4(3) );
% 
% min_x = min([P1(1) P2(1) P3(1) P4(1)]);
% max_x = max([P1(1) P2(1) P3(1) P4(1)]);
% 
% min_y = min([P1(2) P2(2) P3(2) P4(2)]);
% max_y = max([P1(2) P2(2) P3(2) P4(2)]);
% 
% min_z = min([P1(3) P2(3) P3(3) P4(3)]);
% max_z = max([P1(3) P2(3) P3(3) P4(3)]);
% 
% 
% [X,Y] = meshgrid(min_x:max_x, min_y:max_y);
% 
% Z = -(1/normal(3)).*(normal(1).*(X) + normal(2).*(Y) + D);
% 
% % U = [P1;P2;P3;P4];
% % 
% % U_n(1,:) = [X(1,1) Y(1,1) Z(1,1)];
% % U_n(2,:) = [X(1,end) Y(1,end) Z(1,end)];
% % U_n(1,:) = [X(end,1) Y(end,1) Z(end,1)];
% % U_n(1,:) = [X(end,end) Y(end,end) Z(end,end)];
% % 
% % T = maketform('affine3d',U_n,U);
% 
% U = [P1(1:2);P2(1:2);P3(1:2)];
% 
% U_n(1,:) = [X(1,1) Y(1,1)];
% U_n(2,:) = [X(1,end) Y(1,end)];
% U_n(3,:) = [X(end,1) Y(end,1)];
% 
% 
% T = maketform('affine',U_n,U);
% 
% size(T.tdata.T )
% 
% nn = T.tdata.T * [X(end,end) Y(end,end) 1]';
% 
% nn(3) = -(1/normal(3)).*(normal(1).*nn(1) + normal(2).*nn(2) + D);
% 
% figure;
% axis equal
% % subplot(121);imshow(aa,[]);
% % subplot(122);
% cla
% 
% 
% plot3(P1(1),P1(2),P1(3),'r*','MarkerSize',12);hold on
% plot3(P2(1),P2(2),P2(3),'b*','MarkerSize',12);hold on
% plot3(P3(1),P3(2),P3(3),'g*','MarkerSize',12);hold on
% plot3(P4(1),P4(2),P4(3),'k*','MarkerSize',12);hold on
% plot3(nn(1),nn(2),nn(3),'r+','MarkerSize',15);hold on
% 
% surface('XData',X, 'YData',Y, 'ZData',Z, ...
%         'CData',aa, 'CDataMapping','direct', ...
%         'EdgeColor','none', 'FaceColor','texturemap')
    %convert2u8((views.axial(:,:,ax)))
%colormap(gray)

%% Sagittal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = zeros(4,4);

M(1:3,4) = views.sagittal_info{sag}.ImagePositionPatient;
M(4,4) = 1;
M(1:3,1) = views.sagittal_info{sag}.ImageOrientationPatient(1:3).*views.sagittal_info{sag}.PixelSpacing(1);
M(1:3,2) = views.sagittal_info{sag}.ImageOrientationPatient(4:6).*views.sagittal_info{sag}.PixelSpacing(2);



i = [0 511];
j = [0 511];

x = [];
y = [];
z = [];
for k=1:2
    for l=1:2
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
img = ss;  %# Sample image for texture map

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
index = [4 1;3 2];  %# Index used to create 2-by-2 surface coordinates
X1 = x(index);      %# x coordinates of surface
Y1 = y(index);       %# y coordinates of surface
Z1 = z(index);        %# z coordinates of surface

plot3(P1(1),P1(2),P1(3),'r*','MarkerSize',12);hold on
plot3(P2(1),P2(2),P2(3),'b*','MarkerSize',12);hold on
plot3(P3(1),P3(2),P3(3),'g*','MarkerSize',12);hold on
plot3(P4(1),P4(2),P4(3),'k*','MarkerSize',12);hold on

hSurface = surf(X1,Y1,Z1,triTexture,...          %# Plot texture-mapped surface
                'FaceColor','texturemap',...
                'EdgeColor','none');
axis equal  
% normal = cross(P1-P2, P1-P3);
% D = -(normal(1)*P4(1) + normal(2)*P4(2)  + normal(3)*P4(3) );
% 
% % [rows cols] = size(views.coronal(:,:,25));
% % 
% % [X,Z] = meshgrid(P1(3):P2(3), P1(1):P3(1));
% % 
% % Y = -(1/normal(2)).*(normal(3).*Z + normal(1).*X + D);
% 
% min_x = min([P1(1) P2(1) P3(1) P4(1)]);
% max_x = max([P1(1) P2(1) P3(1) P4(1)]);
% 
% min_z = min([P1(3) P2(3) P3(3) P4(3)]);
% max_z = max([P1(3) P2(3) P3(3) P4(3)]);
% 
% [X,Z] = meshgrid(min_x:max_x, min_z:max_z);
% 
% Y = -(1/normal(2)).*(normal(1).*X + normal(3).*Z + D);
% 
% surface('XData',X, 'YData',Y, 'ZData',Z, ...
%         'CData',ss(end:-1:1,end:-1:1,:), 'CDataMapping','direct', ...
%         'EdgeColor','none', 'FaceColor','texturemap')
    %convert2u8(views.sagittal(end:-1:1,end:-1:1,sag))
%colormap(gray)
%% Coronal %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
M = zeros(4,4);

M(1:3,4) = views.coronal_info{cor}.ImagePositionPatient;
M(4,4) = 1;
M(1:3,1) = views.coronal_info{cor}.ImageOrientationPatient(1:3).*views.coronal_info{cor}.PixelSpacing(1);
M(1:3,2) = views.coronal_info{cor}.ImageOrientationPatient(4:6).*views.coronal_info{cor}.PixelSpacing(2);



i = [0 511];
j = [0 511];
x = [];
y = [];
z = [];
for k=1:2
    for l=1:2
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
img = cc;  %# Sample image for texture map

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
%index = [1 4;2 3];  %# Index used to create 2-by-2 surface coordinates
index = [4 1;3 2];
X1 = x(index);      %# x coordinates of surface
Y1 = y(index);       %# y coordinates of surface
Z1 = z(index);        %# z coordinates of surface

plot3(P1(1),P1(2),P1(3),'r*','MarkerSize',12);hold on
plot3(P2(1),P2(2),P2(3),'b*','MarkerSize',12);hold on
plot3(P3(1),P3(2),P3(3),'g*','MarkerSize',12);hold on
plot3(P4(1),P4(2),P4(3),'k*','MarkerSize',12);hold on

hSurface = surf(X1,Y1,Z1,triTexture,...          %# Plot texture-mapped surface
                'FaceColor','texturemap',...
                'EdgeColor','none');
axis equal  
% normal = cross(P1-P2, P1-P3);
% D = -(normal(1)*P4(1) + normal(2)*P4(2)  + normal(3)*P4(3) );
% % 
% % [rows cols] = size(views.coronal(:,:,25));
% % 
% % [Y,Z] = meshgrid(P1(3):P2(3), P1(2):P3(2));
% % 
% % X = -(1/normal(1)).*(normal(3).*Z + normal(2).*Y + D);
% 
% min_y = min([P1(2) P2(2) P3(2) P4(2)]);
% max_y = max([P1(2) P2(2) P3(2) P4(2)]);
% 
% min_z = min([P1(3) P2(3) P3(3) P4(3)]);
% max_z = max([P1(3) P2(3) P3(3) P4(3)]);
% 
% [Y,Z] = meshgrid(min_y:max_y, min_z:max_z);
% 
% X = -(1/normal(1)).*(normal(2).*Y + normal(3).*Z + D);
% 
% surface('XData',X, 'YData',Y, 'ZData',Z, ...
%         'CData',cc, 'CDataMapping','direct', ...
%         'EdgeColor','none', 'FaceColor','texturemap')
%     %convert2u8(views.coronal(end:-1:1,:,cor))
% colormap(gray)
%surf(X, Y, Z, convert2u8(views.coronal(:,:,25)), 'FaceColor', 'texturemap')



%view(3)