warning off
addpath tools

clc
% clear

% load data
% filename='brain.nii';
% V=load_untouch_nii(filename);
% ima=double(V.img);
% ima=ima(:,:,30)';
% res=V.hdr.dime.pixdim(2:4);

ima = double(views.sagittal);
res = [1 1 10];

% filter data
tic
[fima,B]=biascorrector3Dv2(ima,1,4,res);
toc

%plot the results
colormap(gray)
clf
subplot(2,2,1),imagesc(ima);
subplot(2,2,2),imagesc(fima);
subplot(2,2,3),imagesc(B);
h1=histogram(ima);
h2=histogram(fima);
subplot(2,2,4),plot(h1),hold on,plot(h2,'r');
legend('Original','Corrected')


