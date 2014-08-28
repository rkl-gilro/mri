function [m_element,k_element,f_element]=ElementMatrix(e,aantknooppuntenperelement,DegreesOfFreedom,x,y,z,CM,D,External3DElementForce,rho)
% This function computes the element mass, stiffness and the external element volume load.
% Initialize element mass-matrix, and element volume load vector===========
m_element=zeros(aantknooppuntenperelement*DegreesOfFreedom,aantknooppuntenperelement*DegreesOfFreedom); 
f_element=zeros(aantknooppuntenperelement*DegreesOfFreedom,1);                              
% Determine coördinates of nodes of the element============================
ElementCoordinates=[x(CM(:,e)); y(CM(:,e)); z(CM(:,e))]';
%==========================================================================
% Page 292, Eugenio Oñate==================================================
xi_Gausspoint=1/4;
eta_Gausspoint=1/4;
mu_Gausspoint=1/4;
weights=1/6;
% Computation of element matrices==========================================
for i=1:numel(weights)
   xi=xi_Gausspoint(i);            
   eta=eta_Gausspoint(i);
   mu=mu_Gausspoint(i);       
   %=======================================================================
   N=ShapeFunctions(xi,eta,mu);                  
   [B,detJ]=StrainMatrix(ElementCoordinates);
   %ElementMassMatrix======================================================
   m_element=m_element+rho*(N')*N*detJ*weights(i);            
   %Element Body Force=====================================================  
   f_element=f_element+N'*External3DElementForce(:,e)*detJ*weights(i);
end
% f_element yields 1/4*[f3D_x f3D_y f3D_z]*volume_element, page 263, Eugenio Oñate.
% Due to the asymmetrical mesh using tetrahedrons, the global volume loadvector isn't
% symmetric. Applying only volume forces in the z-direction can cause deflections in x and y directions.
% ElementStifnessMatrix====================================================
k_element=(B')*D*B*1/6*detJ;  
% Clarification:(B')*D*B*V_e (detJ = 6V_e)=================================






