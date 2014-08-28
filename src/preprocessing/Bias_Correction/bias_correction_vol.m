function [out, B] = bias_correction_vol(ima, res, umbral)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Bias field correction of a given volume 
%%
%% Inputs:  1. ima -> volume to be filtered
%%          2. res -> voxel resolution (e.g. [1,1,3])
%%          3. umbral -> Background elimination (0:Dont erase, 1: Erase)
%% Outputs: 1. fima -> filtered volume
%%          2. B -> estimated bias field
%%        
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


[out, B] = biascorrector3Dv2(ima, umbral, 4, res);
