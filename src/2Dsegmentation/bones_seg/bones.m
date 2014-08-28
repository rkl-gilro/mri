%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Bones Segmentation First approach:
%%  
%%  Work only with axial view, focus on the first slice
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc
close all
%% Segment the bones:
%%    1. use threshold
%%    2. morphological operations
%%    3. initial 
%%

%% Read the image and denoise it
im = views.axial(:,:,13);
im_den = anisodiff2D(im,10,1/7,30,2);

%% Gradient map


%% Select the ROI
figure;
imshow(im_den,[]);


%% Thresholding
im_thr = zeros(size(im));
im_thr(im_den < 490 & im_den > 350) = 1;

figure;
imshow(im_thr,[]);

%% Morphological operations
% Erode & Fill
fill1 = imfill(im_thr,'holes');
erode = imerode(fill1,ones(5));
fill = imfill(erode,'holes');

figure;
imshow(fill,[]);

% # Connected components
CC = bwconncomp(fill);
numPixels = cellfun(@numel,CC.PixelIdxList);
%[sorted,idx] = sort(numPixels,'ascend');
tmp = fill;

for i = 1:length(numPixels)
    if numPixels(i) < 500
        tmp(CC.PixelIdxList{i}) = 0;
    end
end

figure;
imshow(tmp,[]);

% Skeleton
sk = bwmorph(tmp,'skel',Inf);

or = regionprops(sk,'Orientation','MajorAxisLength', ...
    'MinorAxisLength', 'Eccentricity','Centroid');


[or_label num] = bwlabel(sk);

figure;
imshow(sk,[])
hold on

phi = linspace(0,2*pi,50);
cosphi = cos(phi);
sinphi = sin(phi);

final_comp = [];
final_cen = [];

for k = 1:length(or)
    xbar = or(k).Centroid(1);
    ybar = or(k).Centroid(2);

    a(i) = or(k).MajorAxisLength/2;
    b = or(k).MinorAxisLength/2;

    theta(i) = pi*or(k).Orientation/180;
    
    % Discard orientations
    if ((theta(i) < 3*pi/4 && theta(i) > pi/4) || (theta(i) > -3*pi/4 && theta(i) < -pi/4))
        final_comp = [final_comp theta(i)];
        final_cen = [final_cen;or(k).Centroid];
    end
    
    R = [ cos(theta(i))   sin(theta(i))
         -sin(theta(i))   cos(theta(i))];

    xy = [a(i)*cosphi; b*sinphi];
    xy = R*xy;

    x = xy(1,:) + xbar;
    y = xy(2,:) + ybar;
    
    %polar([0 theta],[0 a],'b+');hold on
    plot(x,y,'r','LineWidth',2);
end
hold off



figure;
imshow(im_den,[]);hold on;
plot(final_cen(:,1),final_cen(:,2),'r*');

% Use Region Growing given the calculated seeds
reg_maxdist = 90;

for i = 1:length(final_cen)
    figure;
    rg = regiongrowing(im_den,round(final_cen(i,2)),round(final_cen(i,1)),reg_maxdist);
    imshow(rg,[]);hold on
end
