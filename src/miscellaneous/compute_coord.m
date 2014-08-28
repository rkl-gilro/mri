function [i1, j1, i2, j2, real_v]  = compute_coord(M1, p1, rows, cols)

real_v = M1 * p1'; % 3D point to 2D point in the frame coordinates

fl  = floor(real_v(1) + 1);
fl2 = floor(real_v(2) + 1);
cl  = ceil(real_v(1) + 1);
cl2 = ceil(real_v(2) + 1);


i1 = min(max(fl2,1),rows);
i2 = min(max(cl2,1),rows);
j1 = min(max(fl,1), cols);
j2 = min(max(cl,1), cols);

