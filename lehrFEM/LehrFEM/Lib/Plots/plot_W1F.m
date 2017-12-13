function varargout = plot_W1F(U,Mesh)
% PLOT_W1F Plot routine for W1F elements.
%
%   FIG = PLOT_W1F(U,MESH) generates a plot of the velocity field
%   represents by the W1F solution U on the struct MESH and returns the 
%   handle FIG to the figure.
%
%   The struct should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh, where
%                M is equal to the number of vertices contained in the
%                mesh.
%    ELEMENTS    M-by-3 matrix specifying the elements of the mesh, where M
%                is equal to the number of elements contained in the mesh.
%    EDGES       P-by-2 matrix specifying the edges of the mesh.
%    VERT2EDGE   M-by-M sparse matrix which specifies whether the two
%                vertices i and j are connected by an edge with number
%                VERT2EDGE(i,j).
%
%   Example:
%
%   fig = plot_W1F(U,Mesh);

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  TOL = .00001;
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  X = Mesh.Coordinates(:,1);
  Y = Mesh.Coordinates(:,2);
  XMin = min(X);
  XMax = max(X);
  YMin = min(Y);
  YMax = max(Y);
  scale = max(XMax-XMin,YMax-YMin);
  meshwidth = get_MeshWidth(Mesh);
  N = max(5,2*floor(scale/meshwidth));
  h1 = scale/(2*N+1);
  h2 = scale/(2*N+3);
  
  % Create container
  
  x = [];
  y = [];
  u = [];
  v = [];
  
  % Calculate value
  
  for i = 1:nElements

      vidx = Mesh.Elements(i,:);
      Vertice = Mesh.Coordinates(vidx,:);
      RU = max(Vertice);
      LD = min(Vertice);
      xMax = RU(1);
      yMax = RU(2);
      xMin = LD(1);
      yMin = LD(2);
      
      P1 = Vertice(1,:);
      P2 = Vertice(2,:);
      P3 = Vertice(3,:);
      bK = P1;
      BK = [P2 - P1;P3 - P1];
      inv_BK = inv(BK);
      TK = transpose(inv_BK);
      det_BK = abs(det(BK));
      
            
      XL = floor((xMin-XMin)/h1);
      XH = ceil((xMax-XMin)/h1);
      YL = floor((yMin-YMin)/h2);
      YH = ceil((yMax-YMin)/h2);
      [X,Y] = meshgrid(XL:XH,YL:YH);
      X = X(:)*h1 + XMin;
      Y = Y(:)*h2 + YMin;
      loc = find((P1(1)-X).*(P2(2)-Y)-(P1(2)-Y).*(P2(1)-X)>TOL &...
                 (P2(1)-X).*(P3(2)-Y)-(P2(2)-Y).*(P3(1)-X)>TOL &...
                 (P3(1)-X).*(P1(2)-Y)-(P3(2)-Y).*(P1(1)-X)>TOL);
      X = X(loc);
      Y = Y(loc);
      ps3 =[X Y];
      
      % Evaluate shape functions at the barycenter of the reference element

      if(~isempty(ps3))

          ps3 = (ps3-ones(size(ps3,1),1)*bK)*inv_BK;
          N = shap_W1F(ps3);
          NS = zeros(size(N));
          NS(:,1:2) = N(:,1:2)*TK;
          NS(:,3:4) = N(:,3:4)*TK;
          NS(:,5:6) = N(:,5:6)*TK;

          % Compute velocity field at barycenters
   
          eidx = [Mesh.Vert2Edge(Mesh.Elements(i,2),Mesh.Elements(i,3)) ...
                  Mesh.Vert2Edge(Mesh.Elements(i,3),Mesh.Elements(i,1)) ...
                  Mesh.Vert2Edge(Mesh.Elements(i,1),Mesh.Elements(i,2))];

          % Determine the orientation

          if(Mesh.Edges(eidx(1),1) == vidx(2))
              p1 = 1;
          else
              p1 = -1;
          end

          if(Mesh.Edges(eidx(2),1) == vidx(3))
              p2 = 1;
          else
              p2 = -1;
          end

          if(Mesh.Edges(eidx(3),1) == vidx(1))
              p3 = 1;
          else
              p3 = -1;
          end

          % Compute velocity field at barycenters

          Lu = U(eidx(1))*p1*NS(:,1) + ...
               U(eidx(2))*p2*NS(:,3) + ...
               U(eidx(3))*p3*NS(:,5);
          Lv = U(eidx(1))*p1*NS(:,2) + ...
               U(eidx(2))*p2*NS(:,4) + ...
               U(eidx(3))*p3*NS(:,6);
          
          u = [u;Lu];
          v = [v;Lv];
          x = [x;X];
          y = [y;Y];

      end
      
  end

  % Generate figure
  
  fig = quiver(x,y,u,v,'b-');
  set(gca,'DataAspectRatio',[1 1 1]);

  if(nargout > 0)
      varargout{1} = fig;
  end
    
return