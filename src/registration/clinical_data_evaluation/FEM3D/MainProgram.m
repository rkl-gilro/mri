% 3D FEM analysis with linear tetrahedral elements. All external
% loads are colored orange. Also a modal analysis of the structure is shown.
clear all;
close all;
hold on;
set(gcf,'color',[1 1 1]);
set(gca,'CameraViewAngle',45,'Projection','perspective');
% Setup mesh===============================================================
x=linspace(0,2,2);
y=linspace(0,8,5);
z=linspace(0,1,2);
[X,Y,Z]=meshgrid(x,y,z);
x=X(:)';
y=Y(:)';
z=Z(:)';
% Connectivity=============================================================
CM=delaunay(X(:),Y(:),Z(:))';
% Set the nodenumbers per element of the six faces=========================
Faces=[1 2 3 1;1 2 4 1;2 3 4 2;3 1 4 3];
% Degrees of freedom of node===============================================
DegreesOfFreedom=3;                                      
% Some factors=============================================================
magnification=8e5; % Displacement magnification
quiver30Dscalefactor=5e-4; % Tame the display-size of the 0D-Boundary forcevectors
quiver31Dscalefactor=1e-2; % Tame the display-size of the 1D-Boundary forcevectors
quiver32Dscalefactor=1e-3; % Tame the display-size of the 2D-Boundary forcevectors
% Define the mode shown in the modal analysis==============================
EigenMode=6;
% Young's modulus, % Poisson's ratio, Constitutive matrix, Density=========
E=2.1e11;                                  
nu=0.3;                                     
D=E/((1+nu)*(1-2*nu))*[1-nu nu   nu   0          0          0          
                       nu   1-nu nu   0          0          0  
                       nu   nu   1-nu 0          0          0
                       0    0    0    (1-2*nu)/2 0          0
                       0    0    0    0          (1-2*nu)/2 0
                       0    0    0    0          0          (1-2*nu)/2];
rho=7.8e3;
% Dofs with predescribed displacement======================================
FixedDofs=zeros(max(CM(:))*DegreesOfFreedom,1);                 % Initialize column with fixed dofs (dofs with prediscribed displacement)
FixedDofsLabeling=zeros(size(FixedDofs));                       % Initialize column with labels with prediscribed displacement
% Label the x,y,z-Dofs of the nodes with prediscribed displacment==========
indicator1=DegreesOfFreedom*[1,6,16,11]-2;                      % x-DOF's with prediscribed displacement
indicator2=DegreesOfFreedom*[1,6,16,11]-1;                      % y-DOF's with prediscribed displacement
indicator3=DegreesOfFreedom*[1,6,16,11];                        % z-DOF's with prediscribed displacement
FixedDofsLabeling([indicator1,indicator2,indicator3])=1;        % Labeling of DOF's with prediscribed displacement
FixedDofs(indicator1)=[0,0,0,0];                                % Displacements of prediscribed x-DOFS 
FixedDofs(indicator2)=[0,0,0,0];                                % Displacements of prediscribed y-DOFS 
FixedDofs(indicator3)=[0,0,0,0];                                % Displacements of prediscribed z-DOFS 
% External 0D boundary forces, like pointforces============================
% [Nodenumber;F_x;F_y;F_z]
External0DElementForce{1}=[10;5000;1000;-1000];
External0DElementForce{2}=[5;0;0;-1000];
% External 1D boundary forces==============================================
% [Node 1 of 1D element;Node 2 of 1D element;F_x;F_y;F_z]
External1DElementForce{1}=[15;20;0;0;50];
% External 2D boundary forces==============================================
% [Node 1 of 2D element;Node 2 of 2D element;Node 3 of 2D element;F_x;F_y;F_z]
External2DElementForce{1}=[13;14;18;0;0;1000];
External2DElementForce{2}=[4;13;14;-2000;0;0];
% External 3D volume loads, like gravity, applies to EACH element==========
% [F_x;F_y;F_z]
External3DElementForce=repmat([0;0;-50],1,size(CM,2));
%==========================================================================
d=zeros(max(CM(:))*DegreesOfFreedom,1);                         % Initialize global displacementvector
f=zeros(numel(d),1);                                            % Initialize global external load vector
M=zeros(numel(d));                                              % Initialize global mass-matrix
K=zeros(numel(d));                                              % Initialize global stiffness-matrix
GlobalNodeStresses=zeros(size(D,1),max(CM(:)));                 % Initialize global stress-matrix
%==========================================================================
CoincidentNodesFrequency=zeros(1,max(CM(:)));
%==========================================================================
% NumberOfFixedDofs stands for the total dofs with prediscribed displacment
NumberOfFixedDofs=numel(indicator1)+numel(indicator2)+numel(indicator3); 
% Keep track of het permutations of all dofs===============================
counter1=0;
counter2=0;
for i=1:numel(d)
   % Check whether a dof is a dof with a prediscribed displacment or not===
   if FixedDofsLabeling(i)==1
      counter1=counter1+1;
      DofPermutationList(i)=counter1; 
      d(counter1)=FixedDofs(i);
    else
      counter2=counter2+1;
      DofPermutationList(i)=NumberOfFixedDofs+counter2; 
    end
end
% Collect the dofs of each element after rearranging with respect to predescribed boundary conditions
for e=1:size(CM,2)
    n=1;
    for j=1:size(CM,1)
        StorageMatrix=DegreesOfFreedom*(CM(j,e)-1);
        for k=1:DegreesOfFreedom
            RearrangedELementDofs(n,e)=DofPermutationList(StorageMatrix+k);
            n=n+1;
        end
    end
end
% Assembly=================================================================
for e=1:size(CM,2)
    [m_element,k_element,fe]=ElementMatrix(e,size(CM,1),DegreesOfFreedom,x,y,z,CM,D,External3DElementForce,rho); 
    Index=RearrangedELementDofs(:,e);
    K(Index,Index)=K(Index,Index)+k_element;
    M(Index,Index)=M(Index,Index)+m_element;
    f(Index)=f(Index)+fe;
end
% Fill up the global external load vector==================================
f=ExternalForces0D1D2D(f,DofPermutationList,External0DElementForce,External1DElementForce,External2DElementForce,x,y,z,DegreesOfFreedom);
% Solve linear system======================================================
d=[d(1:NumberOfFixedDofs)             
     K(NumberOfFixedDofs+1:numel(d),NumberOfFixedDofs+1:numel(d))\...
     (f(NumberOfFixedDofs+1:numel(d))-K(1:NumberOfFixedDofs,NumberOfFixedDofs+1:numel(d))'*...
     d(1:NumberOfFixedDofs))];
% Determine new node coördinates===========================================
Displacement=d(DofPermutationList)*magnification; 
j=1;
for i=1:DegreesOfFreedom:numel(d) 
  x_new(j)=x(j)+Displacement(i);
  y_new(j)=y(j)+Displacement(i+1);
  z_new(j)=z(j)+Displacement(i+2);
  j=j+1;
end
% Compute stresses in the nodes============================================
for e=1:size(CM,2)
   [GlobalNodeStresses]=...
    StressRoutine(d,e,CM,RearrangedELementDofs,x,y,z,D,GlobalNodeStresses);
end
% Plot the undeformed mesh=================================================
for e=1:size(CM,2)
   ElementNodeNumbers=CM(:,e);
   CoincidentNodesFrequency(ElementNodeNumbers)=CoincidentNodesFrequency(ElementNodeNumbers)+ones(1,size(CM,1));
end
Plotmesh(x,y,z,CM,size(CM,2),max(CM(:)),CoincidentNodesFrequency,Faces,'TextmodeOff',[0.75,0.75,0.75]);
% Plot deformed mesh and plot Von Misess Stresses==========================
for e=1:size(CM,2)
   %=======================================================================
   for i=1:size(D,1)
      s{i}=GlobalNodeStresses(i,CM(:,e))./CoincidentNodesFrequency(CM(:,e));
   end
   %=======================================================================          
   VonMisessStress=1./sqrt(2)*sqrt((s{1}-s{2}).^2+(s{2}-s{3}).^2+(s{3}-s{1}).^2+...
       6*s{4}.^2+6*s{5}.^2+6*s{6}.^2); 
   %=======================================================================
   for i=1:size(Faces,1)
      Xnew{i}=x_new(CM(Faces(i,:),e));
      Ynew{i}=y_new(CM(Faces(i,:),e));
      Znew{i}=z_new(CM(Faces(i,:),e));
      ElementNodesStresses{i}=VonMisessStress(Faces(i,:));
      fill3(Xnew{i},Ynew{i},Znew{i},ElementNodesStresses{i});      
   end
end
Plotmesh(x_new,y_new,z_new,CM,size(CM,2),max(CM(:)),CoincidentNodesFrequency,Faces,'TextmodeOn',[0,0,0]);
% Plot 0D boundary forces==================================================
for i=1:numel(External0DElementForce)
   quiver3(x_new(External0DElementForce{i}(1)),y_new(External0DElementForce{i}(1)),z_new(External0DElementForce{i}(1)),...
   External0DElementForce{i}(2),External0DElementForce{i}(3),External0DElementForce{i}(4),...
   quiver30Dscalefactor,'linewidth',2,'Color',[1 0.5 0],'MaxHeadSize',1);
end
% Plot 1D boundary forces==================================================
for i=1:numel(External1DElementForce)
   plot3([x_new(External1DElementForce{i}(1)) x_new(External1DElementForce{i}(2))],...
   [y_new(External1DElementForce{i}(1)) y_new(External1DElementForce{i}(2))],...
   [z_new(External1DElementForce{i}(1)) z_new(External1DElementForce{i}(2))],'Color',[1 0.5 0],'LineWidth',3);
   quiver3(1/2*sum(x_new(External1DElementForce{i}(1:2))),1/2*sum(y_new(External1DElementForce{i}(1:2))),1/2*sum(z_new(External1DElementForce{i}(1:2))),...
   External1DElementForce{i}(3),External1DElementForce{i}(4),External1DElementForce{i}(5),...
   quiver31Dscalefactor,'linewidth',2,'Color',[1 0.5 0],'MaxHeadSize',1);
end
% Plot 2D boundary forces==================================================
for i=1:numel(External2DElementForce)
   fill3(x_new(External2DElementForce{i}(1:3)),y_new(External2DElementForce{i}(1:3)),z_new(External2DElementForce{i}(1:3)),[1 0.5 0]);   
   quiver3(1/3*sum(x_new(External2DElementForce{i}(1:3))),1/3*sum(y_new(External2DElementForce{i}(1:3))),1/3*sum(z_new(External2DElementForce{i}(1:3))),...
   External2DElementForce{i}(4),External2DElementForce{i}(5),External2DElementForce{i}(6),...
   quiver32Dscalefactor,'linewidth',2,'Color',[1 0.5 0],'MaxHeadSize',1);
end
% Plot predescribed boundaries=============================================
for i=1:length(indicator1)
   nodes=(indicator1+2)/DegreesOfFreedom;
   quiver3(x(nodes(i)),y(nodes(i)),z(nodes(i)),1,0,0,'linewidth',1,'Color','m','MaxHeadSize',0.5);  
end;
for i=1:length(indicator2)
   nodes=(indicator2+1)/DegreesOfFreedom;
   quiver3(x(nodes(i)),y(nodes(i)),z(nodes(i)),0,1,0,'linewidth',1,'Color','m','MaxHeadSize',0.5);  
end;
for i=1:length(indicator3)
   nodes=(indicator3)/DegreesOfFreedom;
   quiver3(x(nodes(i)),y(nodes(i)),z(nodes(i)),0,0,1,'linewidth',1,'Color','m','MaxHeadSize',0.5);  
end;
%==========================================================================
light;light;
camlight('left');
lighting gouraud;
title('Von Mises \sigma','FontSize',12','FontName','Times');
box on;
view(100,30);
colorbar;
rotate3d on;
axis equal;
axis off;
camzoom(5);
% Modal analysis===========================================================
M_cut=M(NumberOfFixedDofs+1:end,NumberOfFixedDofs+1:end);
K_cut=K(NumberOfFixedDofs+1:end,NumberOfFixedDofs+1:end);
% Computation of eigenvalues and eigenvectors==============================
[mode,lambda]=eig(M_cut\K_cut,'nobalance');
omega=sqrt(diag(lambda));
[omega,I]=sort(omega,'ascend');
mode=mode(:,I);
%==========================================================================
NodesCoordinates=zeros(numel(d),1);
NodesCoordinates(1:3:end)=x;
NodesCoordinates(2:3:end)=y;
NodesCoordinates(3:3:end)=z;
% Plot eigenmodes==========================================================
mode_coor=zeros(size(K,1),1);
mode_coor(NumberOfFixedDofs+1:end)=mode(:,EigenMode);
%==========================================================================
figure(2)
set(gca,'color',[1 1 1]);
set(gca,'CameraViewAngle',45,'Projection','perspective');
Plotmesh(x_new,y_new,z_new,CM,size(CM,2),max(CM(:)),CoincidentNodesFrequency,Faces,'TextmodeOn',[0,0,0]);
hold on;
frames=20;
MaxStress=0;
for l=1:frames+1,
   cla;  
   T=(l-1)/frames*2*pi;
   node_mode=NodesCoordinates+sin(T)*mode_coor(DofPermutationList);
   % Mode coördinates
   j=1;
   for k=1:DegreesOfFreedom:numel(d)
      x_new(j)=x(j)+node_mode(k);
      y_new(j)=y(j)+node_mode(k+1);
      z_new(j)=z(j)+node_mode(k+2);
      j=j+1;
   end
   %=======================================================================
   GlobalNodeStresses=zeros(size(GlobalNodeStresses));
   %=======================================================================
   for e=1:size(CM,2)
      [GlobalNodeStresses]=...
      StressRoutine(sin(T)*mode_coor,e,CM,RearrangedELementDofs,x,y,z,...
      D,GlobalNodeStresses);
   end
   %Plot the mode==========================================================
   for e=1:size(CM,2)
      %====================================================================
      for i=1:size(D,1)
         s{i}= GlobalNodeStresses(i,CM(:,e))./CoincidentNodesFrequency(CM(:,e));
      end
      %====================================================================             
      VonMisessStress=1./sqrt(2)*sqrt((s{1}-s{2}).^2+(s{2}-s{3}).^2+(s{3}-s{1}).^2+...
          6*s{4}.^2+6*s{5}.^2+6*s{6}.^2);  
      %====================================================================
      for i=1:size(Faces,1)
         Xnew{i}=x_new(CM(Faces(i,:),e));
         Ynew{i}=y_new(CM(Faces(i,:),e));
         Znew{i}=z_new(CM(Faces(i,:),e));
         ElementNodesStresses{i}=VonMisessStress(Faces(i,:));
         fill3(Xnew{i},Ynew{i},Znew{i},ElementNodesStresses{i});
      end
   end
   if l==1
      axis equal;
      axis off;
      view(100,30);
      axis vis3d;
      camzoom(8);
      Ax=axis;
      % The coloring range should be adjusted to different modes, in order to get nice gradients
      caxis([0 2e10]);
      drawnow;
   end
   axis(Ax);
   drawnow;
end



