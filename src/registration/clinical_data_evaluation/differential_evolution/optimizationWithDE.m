function [xmin, tfxmin] = optimizationWithDE(show_output, n_points, ...
    target,outer_shift_cart_op,movements_min, movements_max, movements_mean, PENALTIES, generations, populationSize, name_fun)

%% DIFFERENTIAL EVOLUTION PARAMETERS
S_struct.I_NP         = 50; % 64 150
S_struct.F_weight     = 0.7; % 0.7
S_struct.F_CR         = 0.9;%0.8 o 0.9?? % change
S_struct.I_D          = n_points; % dimensionalidad ;
S_struct.FVr_minbound = -20;
S_struct.FVr_maxbound =  20;
S_struct.I_bnd_constr = 1;
S_struct.I_itermax    = 100; %200 400
S_struct.F_VTR        = 0; % valor muy pequenito que, de alcanzarse, se para el algoritmo
S_struct.I_strategy   = 2; %2 = local-to-best
if show_output > 0
    S_struct.I_refresh    = 10;
else
    S_struct.I_refresh    = 0;
end
S_struct.I_plotting   = 1;

%de_find(0,target,outer_shift_cart_op,movements_min, movements_max, movements_mean, PENALTIES);%[27 0.055 1]);
[xmin, tfxmin] = deopt(name_fun, S_struct); % deopt_find_mask myfun_unc2