function [Pi0 Pi1 Pi2]=BConV(N)
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
%Input: N: the size of IBS.
%Output: Pi0, Pi1, Pi2 to be used in RegMatrix3D. 

%It computes the convolutions between two B-Spline basis functions/ their
%derrivatives. These terms is used to construct the regularization matrix.

B=[-1 3 -3 1; 3 -6 3 0; -3 0 3 0; 1 4 1 0]/6;
%[t^3 t^2 t 1] will be multiplied from the left hand side!
%so its derivaive will be as follows:
dB=[0 0 0 0; -3 9 -9 3; 6 -12 6 0; -3 0 3 0]/6;
d2B=[0 0 0 0; 0 0 0 0; -6 18 -18 6; 6 -12 6 0]/6;
%constructing 4x4 matrices made out of the blending functions.
pi0=INTofMULT(B,B);
pi1=INTofMULT(dB,dB);
pi2=INTofMULT(d2B,d2B);

Pi0=zeros(N); Pi1=zeros(N); Pi2=zeros(N);
%We start moving on the active region (the interval [0 1]).
stepX=1/(N-3);
for i=1:N-3
    %in each subsection there will be 4 active basis in one side
    %to be convolved with 4 other in the other side.
    for r=0:3        
        for s=0:3
            %b_r is related to i+r & b_s is related to i+s.
            Pi0(i+r,i+s)=Pi0(i+r,i+s)+pi0(r+1,s+1)*stepX;
            Pi1(i+r,i+s)=Pi1(i+r,i+s)+pi1(r+1,s+1)*stepX;
            Pi2(i+r,i+s)=Pi2(i+r,i+s)+pi2(r+1,s+1)*stepX;
        end
    end
end
 

function C=INTofMULT(A,B)
%%
C=zeros(4);
for i=1:4
    for j=1:4
        %consider two different columns (showing two basised):
        W=MULT(A(:,i),B(:,j));
        sum=0;
        for k=1:7
            sum=sum+W(k)/(8-k); %it used to be /(9-k) by mistake!
        end
        C(i,j)=sum;
    end
end

function W=MULT(U,V)
%%
%algebraic product of two columnvectors:

W=zeros(7,1);
for i=1:4
    for j=1:4
        %consider the term U(i) & V(j) which is coef. of t^(4-i) & t^(4-j)
        %so the multiply has the power (8-i-j) and must be saved in C(i+j-1)
        W(i+j-1)=W(i+j-1)+U(i)*V(j);
    end
end
