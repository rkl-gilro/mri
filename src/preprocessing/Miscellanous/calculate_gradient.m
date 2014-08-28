function [out_mag out_tan] = calculate_gradient(im)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%%  Calculates the magnitude and gradient direction of a given image
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
[Ix Iy] = gradient(im);

out_mag = sqrt(Ix.^2 + Iy.^2);
out_tan = atan2(Iy,Ix);

figure;
subplot(121);imshow(out_mag,[]);
subplot(122);imshow(out_tan,[]);