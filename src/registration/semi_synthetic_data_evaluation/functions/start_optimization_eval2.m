%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Same as start_optimization_eval
%% In this case the boundary of the grid does not affect the minimization
%% problem, it remains invariant (does not move)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
global source_tri_v
global target_tri_ax
global target_tri_sag
global target_tri_cor

global optimizer
global scalar_v
global scalar_v2
global t_reg

global out_ind

scalar_v = [];
scalar_v2 = [];

preparing_eval;

tic
%% initialization
% mesh0 = [source_tri_v.X(:,1:2);source_tri_v.X(:,2:3);source_tri_v.X(:,1) source_tri_v.X(:,3)]; % [source_tri_v.X(:,1:2);source_tri_v.X(:,2:3)]

tmp1  = source_tri_v.X(out_ind==0,1:2)';
tmp12 = tmp1(:);
tmp2  = source_tri_v.X(out_ind==0,2:3)';
tmp22 = tmp2(:);
tmp3  = [source_tri_v.X(out_ind==0,1) source_tri_v.X(out_ind==0,3)]';
tmp32 = tmp3(:);
X0 = [tmp12' tmp22' tmp32'];
        
options = optimset('LargeScale','off','GradObj','on','TolFun', optimizer.tolfun,'TolX', optimizer.tolx,'Display','iter','MaxIter', optimizer.maxiter,'OutputFcn', {@outfun,@outfun2});
% options = optimset('Display','iter','MaxIter', 100, 'TolFun', .1); %, 'PlotFcns',@optimplotfval);
[xfinal fval exitflag output] = fminunc(@myfun_unc_ortho_eval2, X0', options);

mesh1 = xfinal(1:length(xfinal)/3);
mesh2 = xfinal(length(xfinal)/3 + 1 :2 * length(xfinal)/3 );
mesh3 = xfinal(2*length(xfinal)/3 + 1 :end );

mesh11 = reshape(mesh1',2,length(xfinal)/6);
mesh21 = reshape(mesh2',2,length(xfinal)/6);
mesh31 = reshape(mesh3',2,length(xfinal)/6);

%% Define a new mesh of points first intersection
tmp_points = zeros(size(source_tri_v.X));
tmp_points(out_ind==1,:) = source_tri_v.X(out_ind==1,:);
tmp_points(out_ind==0,:) = [mesh11' source_tri_v.X(out_ind==0,3)]; 

target_tri_ax = TriRep( source_tri_v.Triangulation,[tmp_points(:,1) tmp_points(:,2) tmp_points(:,3)]); % define the new mesh for the axial

%% Define a new mesh of points second intersection
tmp_points = zeros(size(source_tri_v.X));
tmp_points(out_ind==1,:) = source_tri_v.X(out_ind==1,:);
tmp_points(out_ind==0,:) = [source_tri_v.X(out_ind==0,1) mesh21']; 

target_tri_sag = TriRep( source_tri_v.Triangulation,[tmp_points(:,1) tmp_points(:,2) tmp_points(:,3)]); % define the new mesh for the sagittal

%% Define a new mesh of points third intersection
tmp_points = zeros(size(source_tri_v.X));
tmp_points(out_ind==1,:) = source_tri_v.X(out_ind==1,:);
tmp_points(out_ind==0,:) = [mesh31(1,:)' source_tri_v.X(out_ind==0,2) mesh31(2,:)']; 

target_tri_cor = TriRep( source_tri_v.Triangulation,[tmp_points(:,1) tmp_points(:,2) tmp_points(:,3)]); % define the new mesh for the coronal

%% Define the triangulations
% target_tri_ax  = TriRep(source_tri_v.Triangulation, [mesh11' source_tri_v.X(:,3)]);
% target_tri_sag = TriRep(source_tri_v.Triangulation, [source_tri_v.X(:,1) mesh21']);
% target_tri_cor = TriRep(source_tri_v.Triangulation, [mesh31(1,:)' source_tri_v.X(:,2) mesh31(2,:)']);

t_reg = toc