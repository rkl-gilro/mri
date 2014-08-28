function IBSLevelSurf(P,ch,step)
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
%Input: P: IBS coefficients;  
%       ch: surface color; step: the discretization step.
%Output: the IBS surface

run_mine=0
boundx=1; boundy=1; boundz=1;
%---------------------------------------------------------------
%[Dfc,Dfx,Dfy,Dfxx,Dfxy,Dfyx,Dfyy]=giveMoreComponent(c);
%Dfc=@(x,y,z) [x^4 x^3*y x^2*y^2 x*y^3 y^4 x^3 x^2*y x*y^2 y^3 x^2 x*y y^2 x y 1];
%---------------------------------------------------------------
%standard format:
%c=[zeros(1,22),1,.4,.5, 1 0 0 1 0 1 0 0 0 -.5]; c(20)=.5; c(21)=-.7;
%step=.03;
%step=.05;
%figure(1); hold off;
%surface plot:---------
%Xrange=-0:.1:10; Yrange=-0:.1:10;
%box=[.2 .8 .2 .8 .2 .8]
% box=[.1 .9 .1 .9 .1 .9]
box=[0 1 0 1 0 1];
%box=[ .1 .85 .4 .6 .15 .85]
%Xrange=0:step:boundx; Yrange=0:step:boundy; Zrange=0:step:boundz;
Xrange=box(1):step:box(2); Yrange=box(3):step:box(4); Zrange=box(5):step:box(6);

%Xrange=-0:step:boundx; Yrange=-0-.2:step:boundy; Zrange=-0:step:boundz;
%Xrange=-.8:step:.8; Yrange=-.6:step:1; Zrange=.45:step:1;
[X,Y,Z] = meshgrid(Xrange,Yrange,Zrange);
hold on;
w1=waitbar(0,'The mesh is being constructed');
for i=1:size(X,1)%length(Xrange)-1
    waitbar(i/size(X,1),w1);
    for j=1:size(X,2)%length(Yrange)-1
        for k=1:size(X,3)%length(Zrange)-1
            W(i,j,k)=BSpline3D(P,X(i,j,k),Y(i,j,k),Z(i,j,k)); %BSpline3D is another alternative
%             if W(i,j,k)>0
%                 scatter3(X(i,j,k),Y(i,j,k),Z(i,j,k),'r.')
%             else
%                 scatter3(X(i,j,k),Y(i,j,k),Z(i,j,k),'b.')
%             end
        end
    end
end
close(w1);

%Isosurface from MatLab:--------------------
%W=postprocess(W,X,Y,Z,x,y,z,5);
level=0.0001;
p = patch(isosurface(X, Y, Z, W, level));
isonormals(X,Y,Z,W, p)
set(p, 'FaceColor', ch, 'EdgeColor', 'none');
daspect([1 1 1])
view(3)
camlight; lighting phong
%My own function without rendering:--------------------
if run_mine
    figure
    cntr=0;
    for i=2:size(W,1)-1
        levelNo(i)=cntr;
        if cntr==0
            firstInd=i;
        end
        for j=2:size(W,2)-1
            for k=2:size(W,3)-1
                if (W(i,j,k)*W(i-1,j,k)<0)||(W(i,j,k)*W(i,j-1,k)<0)||(W(i,j,k)*W(i,j,k-1)<0)
                    if (W(i,j,k-1)*W(i,j,k+1)<0)
                        cntr=cntr+1;
                        %scatter3(X(i,j,k),Y(i,j,k),Z(i,j,k),'r.')
                        x(cntr)=X(i,j,k); y(cntr)=Y(i,j,k); z(cntr)=Z(i,j,k);
                        %myPatch(c,step,x(cntr),y(cntr),z(cntr));
                    end
                end
            end
        end
    end
end
%box on;
axis(box);
%axis([0,1,0,1,0,1]);
set(gca,'xtick',[ ])
set(gca,'ytick',[ ])
set(gca,'ztick',[ ])
axis off; axis equal;

%removing zeros from the head:
function myTriangle(normal,P1,P2,P3)
hold on;
% normal=cross(P1-P2,P1-P3);
% if norm(normal)
%     normal=abs(normal)/norm(normal);
% else
%     normal=[0 0 0];
% end
col=abs(normal*[1 0 0]')*[1 0 0]
%normal*[1 0 0]'*[1 0 0]
%col=col/norm(col);
if isnan(col)
    col=[1 0 0];
end
X=[P1(1),P2(1),P3(1)]; Y=[P1(2),P2(2),P3(2)]; Z=[P1(3),P2(3),P3(3)];
patch(X,Y,Z,abs(col));



function myPatch(c,step,xi,yi,zi)
%first of all, we should find 4 points on the corners:
[f, fx,fy,fz,fxx,fxy,fxz,fyy,fyz,fzz]=Poly(deg,c,xi,yi,zi);
%I hope fz is nonzero!
DeltaZ1=(-fx/fz-fy/fz)*step/2;
DeltaZ2=(+fx/fz-fy/fz)*step/2;
swapxz=0;
if abs(DeltaZ1)+abs(DeltaZ2)>8*step
    DeltaX1=(-fz/fx-fy/fx)*step/2;
    DeltaX2=(+fz/fx-fy/fx)*step/2;
    %swap x & z:
    swapxz=1;
    tt=xi; xi=zi; zi=tt;
    tt=fx; fx=fz; fz=tt;
    DeltaZ1=DeltaX1; DeltaZ2=DeltaX2;    
end
P=[xi,yi,zi];
normal=[fx, fy, fz]; normal=normal/norm(normal);
%quiver3(P(1),P(2),P(3),n(1),n(2),n(3),'g','LineWidth',2);
P1=[xi+step/2, yi+step/2,zi+DeltaZ1];
P2=[xi-step/2, yi+step/2,zi+DeltaZ2];
P3=[xi-step/2, yi-step/2,zi-DeltaZ1];
P4=[xi+step/2, yi-step/2,zi-DeltaZ2];
if swapxz
    tt=P(1); P(1)=P(3); P(3)=tt;
    tt=P1(1); P1(1)=P1(3); P1(3)=tt;
    tt=P2(1); P2(1)=P2(3); P2(3)=tt;
    tt=P3(1); P3(1)=P3(3); P3(3)=tt;
    tt=P4(1); P4(1)=P4(3); P4(3)=tt;
end
%the extreme case:
%draw all 4 triangles:
myTriangle(normal,P,P1,P2);
myTriangle(normal,P,P2,P3);
myTriangle(normal,P,P3,P4);
myTriangle(normal,P,P4,P1);
% scatter3(P(1),P(2),P(3),120,'filled')
% scatter3(P1(1),P1(2),P1(3),'filled')
% scatter3(P2(1),P2(2),P2(3),'filled')
% scatter3(P3(1),P3(2),P3(3),'filled')
% scatter3(P4(1),P4(2),P4(3),'filled')