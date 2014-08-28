function H=RegMatrix3D(N)
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
%Output: H: IBS Regularization matrix (p'Hp: regularization term).

%This function makes the regularization matrix needed for the L3.
%The output R is a matrix showing the coefficients to minimze the bending
%energy. 
%we start moving on the whole net, and making the term for the
%corresponding control parameters.
%H=zeros(N^3,N^3);
H=sparse(N^3,N^3);
[Pi0 Pi1 Pi2]=BConV(N);
w1=waitbar(0,'Constructing the regularization matrix');
for ik=1:N
    waitbar(ik/N,w1);  
    for jk=1:N
        for kk=1:N
            %-------            
            %first index: row
            k=(ik-1)*N^2+(jk-1)*N+kk;
            for il=1:N
                for jl=1:N
                    for kl=1:N
                        %second index: column
                        l=(il-1)*N^2+(jl-1)*N+kl;                
                        %Here it is the regularization matrix:
                        H(k,l)=Pi2(ik,il)*Pi0(jk,jl)*Pi0(kk,kl)...
                              +2*Pi1(ik,il)*Pi1(jk,jl)*Pi0(kk,kl)...
                              +Pi0(ik,il)*Pi2(jk,jl)*Pi0(kk,kl)...
                              +2*Pi1(ik,il)*Pi0(jk,jl)*Pi1(kk,kl)...
                              +2*Pi0(ik,il)*Pi1(jk,jl)*Pi1(kk,kl)...
                              +Pi0(ik,il)*Pi0(jk,jl)*Pi2(kk,kl);
                    end
                end
            end            
            %-------
        end
    end
end
close(w1);