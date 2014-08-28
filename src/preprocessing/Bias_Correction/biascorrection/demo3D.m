warning off
addpath tools

clc
clear

% load data
filename='brain.nii';
V=load_untouch_nii(filename);
ima=double(V.img);
res=V.hdr.dime.pixdim(2:4);

% filter data
tic
[fima,B]=biascorrector3Dv2(ima,1,4,res);
toc

% save filtered data
nfilename=['m_',filename];
V.img=fima;
save_untouch_nii(V, nfilename);

% save field
nfilename=['bias_',filename];
V.img=B;
save_untouch_nii(V, nfilename);

