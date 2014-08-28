function [fima]=ascm(ima,fimau,fimao,h)

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

s = size(fimau);

p(1) = 2^(ceil(log2(s(1))));
p(2) = 2^(ceil(log2(s(2))));
p(3) = 2^(ceil(log2(s(3))));

pad1 = zeros(p(1),p(2),p(3));
pad2=pad1;
pad3=pad1;
pad1(1:s(1),1:s(2),1:s(3)) = fimau(:,:,:);
pad2(1:s(1),1:s(2),1:s(3)) = fimao(:,:,:); 
pad3(1:s(1),1:s(2),1:s(3)) = ima(:,:,:); 

[af, sf] = farras;
w1 = dwt3D(pad1,1,af);
w2 = dwt3D(pad2,1,af);
w3 = dwt3D(pad3,1,af);

%BayeSkrink for Coeff mixing

for i=1:1:7
    
    tmp = w3{1}{i};
    tmp = tmp(1:round((s(1)-1)/2),1:round((s(2)-1)/2),1:round((s(3)-1)/2));
    sigY = std(tmp(:));
    sigX = sigY^2 - h*h;
    if (sigX < 0)
        T=max(abs(w3{1}{i}(:)));
    else
        T = (h*h) / sqrt(sigX);
    
    end;
    
    w3{1}{i} = abs(w3{1}{i});
    
    dist = w3{1}{i} - T;
    dist = exp(-0.01*dist);
    dist = 1./(1+dist);
    
    w3{1}{i} =  dist.*w1{1}{i} + (1-dist).*w2{1}{i};

end
  
w3{2} = w1{2};

fima = idwt3D(w3,1,sf);
fima = fima(1:s(1),1:s(2),1:s(3));
