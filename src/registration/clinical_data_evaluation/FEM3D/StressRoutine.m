function[GlobalNodeStresses]=StressRoutine(d,e,CM,PermutationList,x,y,z,D,GlobalNodeStresses)
% This function computes the six stresses in each node of an element (the same for every
% node per element in case of linear tetrahedral elements), and stores them
% into the matrix called 'GlobalNodeStresses'.
% Coördinates of the element===============================================
C=[x(CM(:,e)); y(CM(:,e)); z(CM(:,e))]';        
% Elementnodes displacements===============================================
d_element=d(PermutationList(:,e));   
%==========================================================================
for i=1:size(CM,1)    
      [B,detJ]=StrainMatrix(C);
      ElementStrains(:,i)=B*d_element;             
end
ElementStresses=D*ElementStrains;
GlobalNodeStresses(:,CM(:,e))=GlobalNodeStresses(:,CM(:,e))+ElementStresses;

