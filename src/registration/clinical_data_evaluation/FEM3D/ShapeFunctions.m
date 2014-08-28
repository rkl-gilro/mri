function N=ShapeFunctions(xi,eta,mu)
% Isoparametric Shape Functions============================================
NN=[1-xi-eta-mu
    xi
    eta
    mu];
%==========================================================================     
N=[];
for i=1:numel(NN)
   Nsmall=[NN(i) 0     0
           0     NN(i) 0
           0     0     NN(i)];
   N=[N,Nsmall];
end
