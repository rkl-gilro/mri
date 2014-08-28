function F = montecarlo(x, y, z)

global optimizer
global t_reg

global target_tri_ax
global target_tri_sag
global target_tri_cor

global source_tri

size_s = size(source_tri.X);

%% Define the initial condition
tmp1  = source_tri.X(:,1:2)' + [x y]';%(repmat([x y],size_s(1),1))';
tmp12 = tmp1(:);
tmp2  = source_tri.X(:,2:3)' + [y z]';%(repmat([y z],size_s(1),1))';
tmp22 = tmp2(:);
tmp3  = [source_tri.X(:,1) source_tri.X(:,3)]' + [x z]';%(repmat([x z],size_s(1),1))';
tmp32 = tmp3(:);
X0 = [tmp12' tmp22' tmp32'];
        
[F, ~] = myfun_unc_orthoJ(X0);

%% Define the triangulations
target_tri_ax  = TriRep(source_tri.Triangulation, [tmp1' source_tri.X(:,3)]);
target_tri_sag = TriRep(source_tri.Triangulation, [source_tri.X(:,1) tmp2']);
target_tri_cor = TriRep(source_tri.Triangulation, [tmp3(1,:)' source_tri.X(:,2) tmp3(2,:)']);

