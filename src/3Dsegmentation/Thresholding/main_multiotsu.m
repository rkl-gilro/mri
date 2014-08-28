clc; close all;

D = patient_or.coronal;

%% Compute thresholding using Otsu method
[c1,c2]=matitk('fomt',[2,128],double(D));

c3 = ones(size(c1)).*255 & ~(c1|c2);

c3 = c3.*255;

seg_D = c1.*0 + c2.*.5 + c3.*1; 

Otsu(seg_D);