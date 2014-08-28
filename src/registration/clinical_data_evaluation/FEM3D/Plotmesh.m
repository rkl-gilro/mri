function[]=Plotmesh(x,y,z,CM,NumberOfElements,totaantknooppuntenmesh,NumberOfCoincidents,Faces,Textmode,Colour)
% This function plots the mesh with or without the global nodenumbers. The
% colors of the plotted node correspond to the coincidence frequency.
set(gcf,'color',[1 1 1]);
hold on;
NodeColors=[[0 1 0];[1 0 0];[0 0 1];[1 1 0];[0 0 0];[0 1 1];[1 0.5 0];[1 0 0.5]];
for i=1:NumberOfElements
   for j=1:size(Faces,1);
      X=x(CM(Faces(j,:),i));
      Y=y(CM(Faces(j,:),i));
      Z=z(CM(Faces(j,:),i));
      plot3(X,Y,Z,'--','LineWidth',1,'Color',Colour);
   end
   if strcmp(Textmode,'TextmodeOn')==1
      for j=1:size(CM,1)
         text(x(CM(j,i)),y(CM(j,i)),z(CM(j,i)),...
         num2str(CM(j,i)),'FontSize',12','FontName','Times','FontAngle','Normal',...
         'FontWeight','Normal','EdgeColor','none','Color',[0,0,0],'HorizontalAlignment','left','VerticalAlignment','bottom');
      end
   end
end
if strcmp(Textmode,'TextmodeOn')==1
   for i=1:totaantknooppuntenmesh
      colour=NodeColors(mod(NumberOfCoincidents(i)-1,length(NodeColors))+1,:);
      plot3(x(i),y(i),z(i),'o','Markersize',5,'MarkerFaceColor',colour,'MarkerEdgeColor','k')
   end
end
