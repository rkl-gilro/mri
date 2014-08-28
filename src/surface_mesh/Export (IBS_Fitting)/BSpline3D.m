function [f index Buvw] =BSpline3D(P,x,y,z)
%% ====================================================================
% Author: Mohammad Rouhani, Morpheo Team, INRIA Rhone Alpes, (2013)
% Email: mohammad.rouhani@inria.fr
% Title: convolutions between two B-Spline basis functions
% Place of publication: Grenoble, France
% Available from: URL
% http://www.iis.ee.ic.ac.uk/~rouhani/mycodes/IBS.rar
%====================================================================
% When using this software, PLEASE ACKNOWLEDGE the effort that went 
% into development BY REFERRING THE PAPER:
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
% Rouhani M. and Sappa A.D., Implicit B-spline fitting using the 3L 
% algorithm, IEEE Conference on on Image Processing (ICIP'11), 2011.
%::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::: 
%% ====================================================================
%About this code:
%Input: P: IBS coefficients; [x,y,z]: point coordinates.
%Output: f: IBS value at the given points; Buvw: IBS monomials.

M=size(P,1); N=size(P,2); O=size(P,3);
stepX=1/(M-3); stepY=1/(N-3); stepZ=1/(O-3);

i=floor(x/stepX)+1;
j=floor(y/stepY)+1;
k=floor(z/stepZ)+1;
u=x/stepX-i+1; 
v=y/stepY-j+1;
w=z/stepZ-k+1;

%index=[i j];
f=0; h=0;

B=[-1 3 -3 1; 3 -6 3 0; -3 0 3 0; 1 4 1 0]/6;
bu=[u^3 u^2 u 1]*B; 
bv=[v^3 v^2 v 1]*B;
bw=[w^3 w^2 w 1]*B;

for ii=0:3
    for jj=0:3
        for kk=0:3             
            if ((i+ii>0) & (j+jj>0)& (k+kk>0))& ((i+ii<=M) & (j+jj<=N)& (k+kk<=O))
                h=h+1; index(h,1)=i+ii; index(h,2)=j+jj; index(h,3)=k+kk;
                Buvw(h)=bu(ii+1)*bv(jj+1)*bw(kk+1);
                f=f+P(i+ii,j+jj,k+kk)*Buvw(h);
            end
        end
    end
end
