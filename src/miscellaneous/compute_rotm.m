function R = compute_rotm(theta)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Compute the rotational matrix given the 3 angles
%%
%% Inputs:  1. theta -> 1x3 vector containing the angles thetax,y,z
%% Outputs: 1. R     -> 3x3 rotation matrix
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Rx = [1 0 0;0 cos(theta(1)) sin(theta(1));0 -sin(theta(1)) cos(theta(1))];

Ry = [cos(theta(2)) 0 sin(theta(2));0 1 0;-sin(theta(2)) 0 cos(theta(2))];

Rz = [cos(theta(3)) sin(theta(3)) 0; -sin(theta(3)) cos(theta(3)) 0;0 0 1];

R = Rx * Ry * Rz;