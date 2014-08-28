function f=ExternalForces0D1D2D(f,PermutationList,External0DElementForce,External1DElementForce,External2DElementForce,x,y,z,DegreesOfFreedom)
% Assembly of external loads. This function stores pointloads (0D), integrates line tractions (1D), surface tractions (2D), and
% assembles (stores) them. Element volumeloads are integrated and assembled together with the
% element stiffness matrices in the routine called 'ElementMatrix'.
% 0D Boundary external loads===============================================
for i=1:numel(External0DElementForce) 
   node1=External0DElementForce{i}(1);
   xdof1=DegreesOfFreedom*node1-2; % xdof of node of 0D element============
   f(PermutationList(xdof1))  =f(PermutationList(xdof1))  +External0DElementForce{i}(2);
   f(PermutationList(xdof1+1))=f(PermutationList(xdof1+1))+External0DElementForce{i}(3);  
   f(PermutationList(xdof1+2))=f(PermutationList(xdof1+2))+External0DElementForce{i}(4);
end
% 1D Boundary external loads===============================================
psi_Gausspoints=0;
weights=2;
for i=1:numel(External1DElementForce)
   % Initialize the external 1D element loadvector=========================
   f_element=zeros(2*DegreesOfFreedom,1); 
   % First  node of 1D element=============================================
   node1=External1DElementForce{i}(1);
   % Second node of 1D element=============================================
   node2=External1DElementForce{i}(2);  
   %=======================================================================
   % Determinant of Jacobian for 1D lineair element========================
   LengthOfElement=norm([x(node2)-x(node1),y(node2)-y(node1),z(node2)-z(node1)]);
   detJ=1/2*LengthOfElement;                    
   for j=1:numel(weights)        
       psi=psi_Gausspoints(j);               
       N=1/2*[1-psi    0     0      1+psi 0     0;       
              0        1-psi 0      0     1+psi 0;         
              0        0     1-psi  0     0     1+psi];  
       f_element=f_element+(N')*...
           External1DElementForce{i}(3:5)*weights(j)*detJ;
   end
   % f_element yields 1/2*[f1D_x f1D_y f1D_z]*length_element===============
   xdof1=DegreesOfFreedom*node1-2;     % xdof of left node of 1D element
   xdof2=DegreesOfFreedom*node2-2;     % xdof of right node of 1D element
   f(PermutationList(xdof1))  =f(PermutationList(xdof1))  +f_element(1);  
   f(PermutationList(xdof1+1))=f(PermutationList(xdof1+1))+f_element(2);  
   f(PermutationList(xdof1+2))=f(PermutationList(xdof1+2))+f_element(3); 
   f(PermutationList(xdof2))  =f(PermutationList(xdof2))  +f_element(4);
   f(PermutationList(xdof2+1))=f(PermutationList(xdof2+1))+f_element(5);
   f(PermutationList(xdof2+2))=f(PermutationList(xdof2+2))+f_element(6);     
end
% 2D Boundary external loads===============================================
% Page 200, Eugenio Oñate==================================================
xi_Gausspoints =1/3;
eta_Gausspoints=1/3;
weights=1/2;
for i=1:numel(External2DElementForce) 
   % Initialize the external 1D element loadvector=========================
   f_element=zeros(3*DegreesOfFreedom,1); 
   % First  node of 2D element=============================================
   node1=External2DElementForce{i}(1);
   % Second node of 2D element=============================================
   node2=External2DElementForce{i}(2);
   % Third node of 2D element==============================================
   node3=External2DElementForce{i}(3); 
   % Determinant of Jacobian for 2D lineair element========================
   % Area of element=======================================================
   A=1/2*norm(cross([x(node3)-x(node1),y(node3)-y(node1),z(node3)-z(node1)],...
       [x(node3)-x(node2),y(node3)-y(node2),z(node3)-z(node2)]));
   % Determinant (twice the area)==========================================
   detJ=2*A;
   for j=1:numel(weights) 
      xi=xi_Gausspoints(j);
      eta=eta_Gausspoints(j);
      % Shape functions for 2D triangular element==========================
      N=[1-xi-eta 0        0        xi 0   0  eta 0   0;   
         0        1-xi-eta 0        0  xi  0  0   eta 0;
         0        0        1-xi-eta 0  0   xi 0   0   eta];
      f_element=f_element+(N')*...
      External2DElementForce{i}(4:6)*weights(j)*detJ;
   end
   % f_element yields 1/3*[f2D_x f2D_y f2D_z]*area_element=================
   xdof1=DegreesOfFreedom*node1-2;     % xdof of node 1 of 2D element
   xdof2=DegreesOfFreedom*node2-2;     % xdof of node 2 of 2D element
   xdof3=DegreesOfFreedom*node3-2;     % xdof of node 3 of 2D element
   f(PermutationList(xdof1))  =f(PermutationList(xdof1))  +f_element(1);  
   f(PermutationList(xdof1+1))=f(PermutationList(xdof1+1))+f_element(2);  
   f(PermutationList(xdof1+2))=f(PermutationList(xdof1+2))+f_element(3); 
   f(PermutationList(xdof2))  =f(PermutationList(xdof2))  +f_element(4);
   f(PermutationList(xdof2+1))=f(PermutationList(xdof2+1))+f_element(5);
   f(PermutationList(xdof2+2))=f(PermutationList(xdof2+2))+f_element(6);
   f(PermutationList(xdof3))  =f(PermutationList(xdof3))  +f_element(7);
   f(PermutationList(xdof3+1))=f(PermutationList(xdof3+1))+f_element(8);
   f(PermutationList(xdof3+2))=f(PermutationList(xdof3+2))+f_element(9);   
end



