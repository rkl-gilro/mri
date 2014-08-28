function R = compute_rotation(axis_rot,theta) 

ux = [          0  -axis_rot(3)  axis_rot(2);...
       axis_rot(3)           0  -axis_rot(1);...
      -axis_rot(2) axis_rot(1)            0];
  
  u_u = [            axis_rot(1)^2 axis_rot(1) * axis_rot(2) axis_rot(1) * axis_rot(3);...
         axis_rot(1) * axis_rot(2)           axis_rot(2)^2   axis_rot(2) * axis_rot(3);...
         axis_rot(1) * axis_rot(3) axis_rot(2) * axis_rot(3)            axis_rot(3)^2];
     
R = eye(3) .* cos(theta)  + ux .* sin(theta) + (1 - cos(theta)) .* u_u;

% if direction == 1
% R = [cos(theta)  0 sin(theta) ;
%       0          1          0 ;
%      -sin(theta) 0 cos(theta)];
% end
% if direction == 2
% R = [1  0 0;
%      0  cos(theta) sin(theta)  ;
%      0 -sin(theta) cos(theta)];
% end
% if direction == 3
% R = [cos(theta)  sin(theta) 0;
%      -sin(theta) cos(theta) 0;
%           0            0    1];
%      
% end