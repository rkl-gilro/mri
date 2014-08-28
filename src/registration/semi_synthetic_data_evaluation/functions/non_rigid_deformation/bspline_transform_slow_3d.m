function Tlocal = bspline_transform_slow_3d(O_trans,Spacing,X)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%
%% 
%% Inputs:  1. O_trans -> 
%%          2. Spacing -> 
%%          3. X -> 
%%
%% Outputs: 1. Tlocal ->
%% 
%% Function is written by D.Kroon University of Twente (August 2010)
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Make row vectors of input coordinates
x2=X(:,1); y2=X(:,2); z2=X(:,3);

% This code calculates for every coordinate in X, the indices of all
% b-spline knots which have influence on the transformation value of
% this point
[m,l,k]=ndgrid(0:3,0:3,0:3); m=m(:)'; l=l(:)'; k=k(:)';
ixs=floor(x2/Spacing(1));
iys=floor(y2/Spacing(2));
izs=floor(z2/Spacing(3));
ix=repmat(floor(ixs),[1 64])+repmat(m,[length(x2) 1]); ix=ix(:);
iy=repmat(floor(iys),[1 64])+repmat(l,[length(y2) 1]); iy=iy(:);
iz=repmat(floor(izs),[1 64])+repmat(k,[length(z2) 1]); iz=iz(:);

% Size of the b-spline grid
s=size(O_trans);

% Points outside the bspline grid are set to the upper corner
Check_bound=(ix<0)|(ix>(s(1)-1))|(iy<0)|(iy>(s(2)-1))|(iz<0)|(iz>(s(3)-1));
ix(Check_bound)=1; iy(Check_bound)=1; iz(Check_bound)=1;
Check_bound_inv=double(~Check_bound);

% Look up the b-spline knot values in neighborhood of the points in (x2,y2)
Cx=O_trans(ix+iy*s(1) +iz*s(1)*s(2) +                    1).*Check_bound_inv;
Cy=O_trans(ix+iy*s(1) +iz*s(1)*s(2) + s(1)*s(2)*s(3)   + 1).*Check_bound_inv;
Cz=O_trans(ix+iy*s(1) +iz*s(1)*s(2) + s(1)*s(2)*s(3)*2 + 1).*Check_bound_inv;
Cx=reshape(Cx,[length(x2) 64]);
Cy=reshape(Cy,[length(x2) 64]);
Cz=reshape(Cz,[length(x2) 64]);

% Calculate the b-spline interpolation constants u,v in the center cell
% range between 0 and 1
v  = (x2-ixs*Spacing(1))/Spacing(1);
u  = (y2-iys*Spacing(2))/Spacing(2);
w  = (z2-izs*Spacing(3))/Spacing(3);

% Get the b-spline coefficients in a matrix W, which contains
% the influence of all knots on the points in (x2,y2)
W=bspline_coefficients(v,u,w);

% Calculate the transformation of the points in (x2,y2) by the b-spline grid
Tlocal(:,1)=sum(W.*Cx,2);
Tlocal(:,2)=sum(W.*Cy,2);
Tlocal(:,3)=sum(W.*Cz,2);