function [fima]=hsm(fimau,fimao)

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

%                          Details on Wavelet mixing                         
%***************************************************************************
 %  The hard wavelet subbands mixing is described in:                      *
 %                                                                         *
 %  P. Coupe, S. Prima, P. Hellier, C. Kervrann, C. Barillot.              *
 %  3D Wavelet Sub-Bands Mixing for Image Denoising                        *
 %  International Journal of Biomedical Imaging, 2008                      * 
% ***************************************************************************/

s = size(fimau);

p(1) = 2^(ceil(log2(s(1))));
p(2) = 2^(ceil(log2(s(2))));
p(3) = 2^(ceil(log2(s(3))));

pad1 = zeros(p(1),p(2),p(3));
pad2 = pad1;
pad1(1:s(1),1:s(2),1:s(3)) = fimau(:,:,:);
pad2(1:s(1),1:s(2),1:s(3)) = fimao(:,:,:); 

[af, sf] = farras;
w1 = dwt3D(pad1,1,af);
w2 = dwt3D(pad2,1,af);

  
  w1{1}{3} = (w2{1}{3}+w1{1}{3})/2;
  w1{1}{5} = (w2{1}{5}+w1{1}{5})/2;
  w1{1}{6} = (w2{1}{6}+w1{1}{6})/2;
  w1{1}{7} = (w2{1}{7}+w1{1}{7})/2;

fima = idwt3D(w1,1,sf);
fima = fima(1:s(1),1:s(2),1:s(3));