%% ====================================================================
% Author: Mohammad Rouhani, Morpheo Team, INRIA Rhone Alpes, (2013)
% Email: mohammad.rouhani@inria.fr
% Title: convolutions between two B-Spline basis functions
% Place of publication: Grenoble, France
% Available from: URL
% http://www.iis.ee.ic.ac.uk/~rouhani/mycodes/IBS.rar
%====================================================================
% When using this software, PLEASE ACKNOWLEDGE the effort that went 
% into development BY REFERRING THE PAPER:
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
% Rouhani M. and Sappa A.D., Implicit B-spline fitting using the 3L 
% algorithm, IEEE Conference on on Image Processing (ICIP'11), 2011.
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
%% ==IP FITTING============================================================
%1. Load the data points:
%load('Bunny.mat'); %load('Duck.mat'); %or any other data set
%Note: the data points must be normalized (in unit cube).
model     = st_all(2).vertices;
model_tri = double(st_all(2).faces);

min_x = min(model(:,1));
min_y = min(model(:,2));
min_z = min(model(:,3));
% 
model(:,1) = model(:,1) + min_x*((-1)^(sign(min_x)));
model(:,2) = model(:,2) + min_y*((-1)^(sign(min_y)));
model(:,3) = model(:,3) + min_z*((-1)^(sign(min_z)));

for i=1:size(model,1)
    
    model(i,1) = model(i,1)/max(model(:,1));
    model(i,2) = model(i,2)/max(model(:,2));
    model(i,3) = model(i,3)/max(model(:,3));
    
end


%2. Call the 3L algorithm to compute IP coefficient vector:
L = 30; %regularization parameter; increase it for a coarser surface.
P = IBSL3_3DTRI(.01, 20, L, model, model_tri); %IBS size can be increased to 30.

%3. Represent the IP surface (zero leve set): 
IBSLevelSurf(P,[.9 .9 .8],0.03); %the visulaization step can be decreased;
%if the surface is not completly reconstruced, change "box" parameters.
