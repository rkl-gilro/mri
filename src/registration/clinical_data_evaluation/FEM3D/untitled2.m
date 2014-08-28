% interiorPoints = rand(200,3);      %# Generate 200 3-D points
% DT = DelaunayTri(interiorPoints);  %# Create the tetrahedral mesh
% hullFacets = convexHull(DT);       %# Find the facets of the convex hull
% 
% %# Plot the scattered points:
% subplot(2,2,1);
% scatter3(interiorPoints(:,1),interiorPoints(:,2),interiorPoints(:,3),'.');
% axis equal;
% title('Interior points');
% 
% %# Plot the tetrahedral mesh:
% subplot(2,2,2);
% tetramesh(DT);
% axis equal;
% title('Tetrahedral mesh');
% 
% %# Plot the 3-D convex hull:
% subplot(2,2,3);
% trisurf(hullFacets,DT.X(:,1),DT.X(:,2),DT.X(:,3),'FaceColor','c')
% axis equal;
% title('Convex hull');

% 
%   t = linspace(0.6,5.7,500)';
%   X = 2*[cos(t),sin(t)] + rand(500,2);
%   subplot(221), alphavol(X,inf,1);
%   subplot(222), alphavol(X,1,1);
%   subplot(223), alphavol(X,0.5,1);
%   subplot(224), alphavol(X,0.2,1);

  [x,y,z] = sphere;
  ii = abs(z) < 0.4;
  X = [x(ii),y(ii),z(ii)];
  X = [X; 0.8*X];
  subplot(211), alphavol(X,inf,1);
  subplot(212), alphavol(X,0.5,1);