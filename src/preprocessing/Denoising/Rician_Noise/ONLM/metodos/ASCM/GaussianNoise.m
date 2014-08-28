% Pierrick Coupe - pierrick.coupe@gmail.com                                  
% Jose V. Manjon - jmanjon@fis.upv.es                                        
% Brain Imaging Center, Montreal Neurological Institute.                     
% Mc Gill University                                                         
%                                                                            
% Copyright (C) 2008 Pierrick Coupe and Jose V. Manjon      

%****************************************************************************
%              3D Adaptive Multiresolution Non-Local Means Filter           *
%            P. Coupe a, J. V. Manjon, M. Robles , D. L. Collin             * 
%****************************************************************************

%                          Details on ONLM filter                        */
%***************************************************************************
 %  The ONLM filter is described in:                                       *
 %                                                                         *
 %  P. Coupe, P. Yger, S. Prima, P. Hellier, C. Kervrann, C. Barillot.     *
 %  An Optimized Blockwise Non Local Means Denoising Filter for 3D Magnetic*
 %  Resonance Images. IEEE Transactions on Medical Imaging, 27(4):425-441, * 
 %  Avril 2008                                                             *
% ***************************************************************************/

%                          Details on Wavelet mixing                       */
%***************************************************************************
 %  The hard wavelet subbands mixing is described in:                      *
 %                                                                         *
 %  P. Coupe, S. Prima, P. Hellier, C. Kervrann, C. Barillot.              *
 %  3D Wavelet Sub-Bands Mixing for Image Denoising                        *
 %  International Journal of Biomedical Imaging, 2008                      * 
% ***************************************************************************/

addpath wavelet

warning off;

clear all;
clc;
close all;
colormap(gray);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% do not forget to compile the mex file
% in the matlab prompt: 
% > mex onlm.c
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Opening the brainweb data
% Chose the modality
name ='data/t1_icbm_normal_1mm_pn0_rf0.rawb';
%  name ='data/t2_icbm_normal_1mm_pn0_rf0.rawb';
% name ='data/pd_icbm_normal_1mm_pn0_rf0.rawb';

fid = fopen(name,'r');    
s=[181,217,181];
ima=zeros(s(1:3));
for z=1:s(3),    
  ima(:,:,z) = fread(fid,s(1:2),'uchar');
end;
fclose(fid);

name ='data/phantom_1.0mm_normal_crisp.rawb';
fid = fopen(name,'r');    
s=[181,217,181];
Label=zeros(s(1:3));
for z=1:s(3),    
  Label(:,:,z) = fread(fid,s(1:2),'uchar');
end;
fclose(fid);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Comment this part to filter the entire image.
ima=ima(:,:,70:100);
Label=Label(:,:,70:100);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


s=size(ima);

% minimum and maximum level of noise in pourcentage
minl = 3;
maxl = 21;

% 3 to 21% of noise
for i=minl:2:maxl

i

% Create noisy data with Gaussian noise
level=i*max(ima(:))/100;
rima=ima+level*randn(s);

% The minimum of the image has to be positive 
% due to the tests on mean and var 
mini = min(rima(:));
rima = rima + abs(mini);

%params
M=3;
alpha=1;


% Filtering with Su parameters: small patch
tic,
h=level;
fima1=onlm(rima,M,alpha,h);
fima1=fima1 - abs(mini);
t1(i)=toc;

% Filtering with So parameters: big patch 
tic,
h=level;
fima2=onlm(rima,M,alpha+1,h);
fima2=fima2 - abs(mini);
t2(i)=toc;


% Hard wavelet Coefficient Mixing
tic,
fima3 = hsm(fima1,fima2);
t3(i)=toc;
t3(i)=t3(i)+t2(i)+t1(i);


% Adaptive wavelet coefficient Mixing
tic,
fima4 = ascm(rima,fima1,fima2,level);
t4(i)=toc;
t4(i)=t4(i)+t2(i)+t1(i);


% Coupe blockwise approach with parameters used in TMI 08
% Search area of 11^3 voxels
tic,
h=level;
fima1=onlm(rima,M+2,alpha,h);
fima1=fima1 - abs(mini);
t1(i)=toc;

rima = rima - abs(mini);

% Remove residual error in reconstructions after IDWT
ind=find(fima3<0);
fima3(ind) = 0;
ind=find(fima4<0);
fima4(ind) = 0;

% Creation of the mask (i.e. without background)
ind=find(Label>0);

% Computation of RMSE
error1(i)=sqrt(mean((fima1(ind)-ima(ind)).^2));
error3(i)=sqrt(mean((fima3(ind)-ima(ind)).^2));
error4(i)=sqrt(mean((fima4(ind)-ima(ind)).^2));

% Computation of PSNR
range=255;
psnr1(i)=20*log10(range/error1(i))
psnr3(i)=20*log10(range/error3(i))
psnr4(i)=20*log10(range/error4(i))

% Figure display
slice = round(s(3)/2);

figure;
c_min1 = 0;
c_max1 = max(ima(:));
c_min2 = min(fima1(:)- rima(:));
c_max2 = max(fima1(:)- rima(:));
colormap(gray);

subplot(2,4,1), imagesc(rot90(ima(:,:,slice)),[c_min1 c_max1]);
title('Ground truth');
subplot(2,4,5), imagesc(rot90(rima(:,:,slice)),[c_min1 c_max1] );
tit = sprintf('Noisy image with %d %% of noise',i);
title(tit);

subplot(2,4,2), imagesc(rot90(fima1(:,:,slice)),[c_min1 c_max1]);
tit = sprintf('ONLM %.2f dB',psnr1(i));
title(tit);
tit = sprintf('Time %.2f s',t1(i));
xlabel(tit);
subplot(2,4,6), imagesc(rot90(fima1(:,:,slice) - rima(:,:,slice)),[c_min2 c_max2]);
title('Residual');

subplot(2,4,3), imagesc(rot90(fima3(:,:,slice)),[c_min1 c_max1]);
tit = sprintf('ONLM with HSM %.2f dB',psnr3(i));
title(tit);
tit = sprintf('Time %.2f s',t3(i));
xlabel(tit);
subplot(2,4,7), imagesc(rot90(fima3(:,:,slice) - rima(:,:,slice)),[c_min2 c_max2]);
title('Residual');


subplot(2,4,4), imagesc(rot90(fima4(:,:,slice)),[c_min1 c_max1]);
tit = sprintf('ONLM with ASCM %.2f dB',psnr4(i));
title(tit);
tit = sprintf('Time %.2f s',t4(i));
xlabel(tit);
subplot(2,4,8), imagesc(rot90(fima4(:,:,slice) - rima(:,:,slice)),[c_min2 c_max2]);
title('Residual');

end


figure;
plot(minl:2:maxl,error1(minl:2:maxl),'k')
 hold on
 plot(minl:2:maxl,error3(minl:2:maxl),'r')
 plot(minl:2:maxl,error4(minl:2:maxl),'g')
title('RMSE')
legend('ONLM','ONLM with HSM','ONLM with ASCM');


figure;
plot(minl:2:maxl,psnr1(minl:2:maxl),'k')
hold on
 plot(minl:2:maxl,psnr3(minl:2:maxl),'r')
 plot(minl:2:maxl,psnr4(minl:2:maxl),'g')
title('PSNR')
legend('ONLM','ONLM with HSM','ONLM with ASCM');


figure;
plot(minl:2:maxl,t1(minl:2:maxl),'k')
 hold on
 plot(minl:2:maxl,t3(minl:2:maxl),'r')
 plot(minl:2:maxl,t4(minl:2:maxl),'g')
title('Time')
legend('ONLM','ONLM with HSM','ONLM with ASCM');




