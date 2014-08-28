function P=IBSL3_3DTRI(r,N,L,pnt,tri)
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
%Input: pnt: points; tri: the triangular mesh; 
%       N: IBS size; r: the offset distance;
%Output: P: the IP coefficient vector

%This 3L algorithm is specially configurated for Implicit BSplines (IBS).
P=zeros(N,N,N);
%---------------------------------------------------------------
x=pnt(:,1);y=pnt(:,2);z=pnt(:,3);

disp('The first step: computing the internal & external offset:')
%r=.4;%the offset distance
%%[s,t]=Offset3D2(r,x,y,z); hold on;
%load('BUNNY_OFF_005.mat');
%[s,t]=NewOffset(r,pnt, tri); hold on;
[s,t]=NewOffset(r, pnt, tri); hold on;

if size(s,2)>3
    s=s'; t=t';
end
x(isnan(x)) = 0;
y(isnan(y)) = 0;
z(isnan(z)) = 0;
s(isnan(s)) = 0;
t(isnan(t)) = 0;
%[s,t]=NewOffset(r,pnt, tri); hold on; s=s'; t=t';
hold on
scatter3(s(:,1),s(:,2),s(:,3),'b','filled');
scatter3(pnt(:,1),pnt(:,2),pnt(:,3),'g','filled');
scatter3(t(:,1),t(:,2),t(:,3),'r','filled');
%--------------------------------
%definition of variables:
n=numel(x);
%M=size(P,1); N=size(P,2); O=size(P,3);
b=-0*ones(n,1);%warning: it must be -1
%MM=zeros(n,N^3);
%define the offset distance and the expected value
eps=1.0; 
b=[b; b+eps; b-eps];
tic
%---------------------------------------------------------------
disp('The second step: Computing the control values of IBS.')
tic;
M0=BlockMatrix3D(P,x,y,z);                M0=sparse(M0);
MP=BlockMatrix3D(P,s(:,1),s(:,2),s(:,3)); MP=sparse(MP);
MN=BlockMatrix3D(P,t(:,1),t(:,2),t(:,3)); MN=sparse(MN);
Time1=toc;
MM=[M0; MP; MN];
clear M0; clear MP; clear MN;

%title(cond(MM));
%c=inv(MM'*MM)*(MM'*b);
%vP=inv(MM'*MM+L*eye(N^3))*(MM'*b); 
if (N==20)
    load('H20.mat');
else if (N==30)
        load('H30.mat');
    else
        H=RegMatrix3D(N);
    end
end
%H=RegMatrix3D(N);
H=sparse(H);
%H=eye(N^3);
H=sparse(H); MM=sparse(MM);
MAT=sparse(MM'*MM+ L*H); bb=(MM'*b);
tic;
vP=MAT\bb; %vP=inv(MM'*MM+ L*H)*(MM'*b); 
Time2=toc 
% This is in the vector form; it should be converted to the matrix form.
for l=1:N^3
    ind=l;
    i=floor((ind-1)/(N^2))+1;
    ind=ind-(i-1)*N^2;
    j=floor((ind-1)/N)+1;
    k=ind-(j-1)*N;
    P(i,j,k)=vP(l);
end
%--------------------------------
disp('The third step: Representing the result.')
disp('Please wait! The surface is being renderd.')
close all
%LevelSurf(P,.05)

% for i=1:M
%     for j=1:N
%         scatter3(i/M,j/N,P(i,j),'filled'); %it is NOT precisely located. 
%     end
% end
%----------------------------------puting additional constraint:
%just one: on p_31
% i=9;j=10;k=11;  W=[Dfc(x(i),y(i)); Dfc(x(j),y(j)); Dfc(x(k),y(k)) ];
% lambda=-2*inv(W*inv(M'*M)*W')*W*inv(M'*M)*M'*b;
% c=inv(M'*M)*(M'*b)+.5*inv(M'*M)*W'*lambda;
numIND=N^3;
hold on; %polyplot(c,-20,20)



function MM=BlockMatrix3D(P,x,y,z)
%%
n=numel(x); N=size(P,1); %N=size(P,2); O=size(P,3); 
%MM=zeros(n,N^3);
MM=sparse(n,N^3);
w1=waitbar(0,'Constructing the block matrix');

for i=1:n
    if (mod(i,200)==0)
        waitbar(i/n,w1);    
    end  

    [f index Buvw]=BSpline3D(P,x(i),y(i),z(i));
    hh=numel(Buvw);
    %ii=index(1); jj=index(2); kk=index(3);
    %hold on;    scatter(index(1),index(2))
    for h=1:hh
        ii=index(h,1); jj=index(h,2); kk=index(h,3);
        %I am defining the rule of vectorization:
        new_index=N^2*(ii-1)+N*(jj-1)+kk;
        %it must be modified back based on this rule!
        MM(i,new_index)=Buvw(h);

    end
    %close all
end
1;
close(w1);
