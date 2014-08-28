function masks = Laws_filter_definition()
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% Compute the 3x3 and 5x5 masks for
%% Laws filtering
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%% 3 x 3 masks %%%%%%%%%%%%%%%%%%%%%%%%%%%
L3 = [1 2 1];
E3 = [1 0 -1];
S3 = [1 -2 1];

masks.L3L3 = L3' * L3;
masks.L3E3 = L3' * E3;
masks.L3S3 = L3' * S3;

masks.E3L3 = E3' * L3;
masks.E3E3 = E3' * E3;
masks.E3S3 = E3' * S3;

masks.S3L3 = S3' * L3;
masks.S3E3 = S3' * E3;
masks.S3S3 = S3' * S3;

%% 5 x 5 masks %%%%%%%%%%%%%%%%%%%%%%%%%%%
L5 = [1 4 6 4 1];
E5 = [-1 -2 0 2 1];
S5 = [-1 0 2 0 -1];
W5 = [-1 2 0 -2 1];
R5 = [1 -4 6 -4 1];

masks.L5L5 = L5'*L5;
masks.L5E5 = L5'*E5;
masks.L5S5 = L5'*S5;
masks.L5W5 = L5'*W5;
masks.L5R5 = L5'*R5;

masks.E5L5 = E5'*L5;
masks.E5E5 = E5'*E5;
masks.E5S5 = E5'*S5;
masks.E5W5 = E5'*W5;
masks.E5R5 = E5'*R5;

masks.S5L5 = S5'*L5;
masks.S5E5 = S5'*E5;
masks.S5S5 = S5'*S5;
masks.S5W5 = S5'*W5;
masks.S5R5 = S5'*R5;

masks.W5L5 = W5'*L5;
masks.W5E5 = W5'*E5;
masks.W5S5 = W5'*S5;
masks.W5W5 = W5'*W5;
masks.W5R5 = W5'*R5;

masks.R5L5 = R5'*L5;
masks.R5E5 = R5'*E5;
masks.R5S5 = R5'*S5;
masks.R5W5 = R5'*W5;
masks.R5R5 = R5'*R5;