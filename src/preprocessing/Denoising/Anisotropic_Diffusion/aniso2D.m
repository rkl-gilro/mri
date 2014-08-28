function nois = aniso2D(vol,num_iter,delta,kappa,func)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Compute anisotropic diffusion given a volume
%% 
%% Inputs:  1. vol -> 3d matrix
%%          2. [num_iter,delta,kappa,func] -> anisotropic diffusion parameters
%%     
%% Outputs: 1. nois -> filtered volume
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% Check parameters %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if nargin < 2
  % default parameters
    num_iter = 20;
    delta = 1/7;
    kappa = 30;
    func = 2;
end

if nargin < 3
    delta = 1/7;
    kappa = 30;
    func = 2;
end

if nargin < 4
   kappa = 30;
   func = 2;
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:size(vol,3)
    nois(:,:,i) = anisodiff2D(vol(:,:,i),num_iter,delta,kappa,func);
end