addpath ./metodos
addpath ./metodos/ASCM
addpath ./metodos/ASCM/wavelet

clc
clear
colormap(gray)

% read volume
name ='t1_icbm_normal_1mm_pn0_rf0.rawb';

fid = fopen(name,'r');    
s=[181,217,181];
ima=zeros(s(1:3));
for z=1:s(3),    
  ima(:,:,z) = fread(fid,s(1:2),'uchar');
end;
fclose(fid);
ima=double(ima);

% select a small part to rapidly check the code (comment this line for full volume processing)
ima=ima(101:150,101:150,101:150);


indi=find(ima>10);
s=size(ima);
R=max(ima(:));
sw = [1 1 1]; 

for i=1:2:9
i
% add noise
level=i*max(ima(:))/100   
rima=sqrt((ima+level*randn(s)).^2+(level*randn(s)).^2);

% Method 1: Block-wise NLM3D Coupe 2008 
tic
v=5;
fima0=MBONLM3D(rima,v,1,level,1);
bnlm3d=toc

% Method 2: Wavelet coefficient Mixing
nv=3;
tic
fima11=MBONLM3D(rima,nv,1,level,1);
fima12=MBONLM3D(rima,nv,2,level,1);
fima1 = hsm(fima11,fima12);
WSM=toc

% Method 3: 3D ODCT filtering (proposed)
tic
fima2=cM_ODCT3D(rima,level,1);
odct3d=toc

% Method 4: PRI-NLM3D filtering (proposed)
mv=5;
fima3=cM_RI_NLM3D(rima,mv,1,level,fima2,1); 
proposed=toc

% results
error0(i)=sqrt(mean((ima(indi)-fima0(indi)).^2))
error1(i)=sqrt(mean((ima(indi)-fima1(indi)).^2))
error2(i)=sqrt(mean((ima(indi)-fima2(indi)).^2))
error3(i)=sqrt(mean((ima(indi)-fima3(indi)).^2))

psnr0(i)=20*log10(R/error0(i))
psnr1(i)=20*log10(R/error1(i))
psnr2(i)=20*log10(R/error2(i))
psnr3(i)=20*log10(R/error3(i))

ssim0(i)= ssim_index3d(fima0,ima,sw,indi)
ssim1(i)= ssim_index3d(fima1,ima,sw,indi)
ssim2(i)= ssim_index3d(fima2,ima,sw,indi)
ssim3(i)= ssim_index3d(fima3,ima,sw,indi)

end

mean(error0)
mean(error1)
mean(error2)
mean(error3)

clf
plot(1:2:i,error0(1:2:i),'bx-')
hold on
plot(1:2:i,error1(1:2:i),'kp-')
plot(1:2:i,error2(1:2:i),'mo-')
plot(1:2:i,error3(1:2:i),'r^-')
xlabel('Noise level (%)')
ylabel('RMSE')
legend('Coupe-Block','WSM','ODCT3D','PRI-NLM3D');

figure
clf
plot(1:2:i,ssim0(1:2:i),'bx-')
hold on
plot(1:2:i,ssim1(1:2:i),'kp-')
plot(1:2:i,ssim2(1:2:i),'mo-')
plot(1:2:i,ssim3(1:2:i),'r^-')
xlabel('Noise level (%)')
ylabel('SSIM')
legend('Coupe-Block','WSM','ODCT3D','PRI-NLM3D');

figure
mv=round(max(ima(:)));
map=zeros(mv,3);
map(:,1)=(0:mv-1)/(mv-1);
map(:,2)=map(:,1);
map(:,3)=map(:,1);
colormap(map);
clf
subplot(3,3,1),image(imrotate(ima(:,:,round(s(3)/2)),90)),xlabel('Original')
subplot(3,3,2),imagesc(imrotate(rima(:,:,round(s(3)/2)),90)),xlabel('Noisy')
subplot(3,3,3),image(imrotate(fima0(:,:,round(s(3)/2)),90)),xlabel('Coupe')
subplot(3,3,4),image(imrotate(fima1(:,:,round(s(3)/2)),90)),xlabel('WSM')
subplot(3,3,5),image(imrotate(fima2(:,:,round(s(3)/2)),90)),xlabel('ODCT3D')
subplot(3,3,6),image(imrotate(fima3(:,:,round(s(3)/2)),90)),xlabel('Proposed')

subplot(3,3,7),imagesc(abs(imrotate(rima(:,:,round(s(3)/2))-fima1(:,:,round(s(3)/2)),90)*1)),xlabel('WSM')
subplot(3,3,8),imagesc(abs(imrotate(rima(:,:,round(s(3)/2))-fima2(:,:,round(s(3)/2)),90)*1)),xlabel('ODCT3D')
subplot(3,3,9),imagesc(abs(imrotate(rima(:,:,round(s(3)/2))-fima3(:,:,round(s(3)/2)),90)*1)),xlabel('Proposed')
 
