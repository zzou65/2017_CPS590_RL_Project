function varargout = plot_Stokes(U,Mesh,type)
% PLOT_STOKES Plot the solution to Stokes problem
%
%   PLOT_STOKES(U,MESH,TYPE) generates a plot of the solution U for the 
%   Stokes problem on the struct MESH for the element pair defined by TYPE:
%    MINIP0 MINI-elements for the velocity field and P0-elements for the
%           pressure.
%    P1P0   P1-elements for the velocity field and P0-elements for the
%           pressure.
%    P2P0   P2-elements for the velocity field and P0-elements for the
%           pressure.
%    CRP0   Crouzeix-Raviart-elements for the velocity field and
%           P0-elements for the pressure.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh, where
%                M is equal to the number of vertices contained in the
%                mesh.
%    ELEMENTS    M-by-3 matrix specifying the elements of the mesh, where
%                M is equal to the number of elements contained in the
%                mesh.
%    EDGES       P-by-2 matrix specifying all edges of the mesh.
%    VERT2EDGE   M-by-M sparse matrix which specifies wheter the two
%                vertices i and j are connected by an edge with number
%                VERT2EDGE(i,j).
%   
%   H = PLOT_STOKES(U,MESH,TYPE) also returns the handle to the figure.
%
%   Example:
%   
%   plot_Stokes(U,Mesh,'MINI');

%   Copyright 2005-2005 Patrick Meury & Kah-Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants

  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  nEdges = size(Mesh.Edges,1);
  
  fig = figure('Name','Pressure and veclocity field');
  if(~ishold)
    hold on;
  end
  if(strcmp(upper(type),'MINIP1'))
      
    % Extract components from solution vector  
      
    U1 = U(1:(nCoordinates+nElements));
    U2 = U(nCoordinates+nElements+(1:(nCoordinates+nElements)));
    P = U(2*(nCoordinates+nElements)+(1:nCoordinates));
    
    % Generate figures
    
    plot_p_P1(P,Mesh);
    plot_u_MINI(U1,U2,Mesh);
    
  elseif(strcmp(upper(type),'P2P0'))
      
    % Extract components from solution vector  
      
    U1 = U(1:(nCoordinates+nEdges));
    U2 = U(nCoordinates+nEdges+(1:(nCoordinates+nEdges)));
    P = U(2*(nCoordinates+nEdges)+(1:nElements));
    
    % Generate figures
    
    plot_p_P0(P,Mesh);  
    plot_u_P2(U1,U2,Mesh);
      
  elseif(strcmp(upper(type),'P1P0'))
      
    % Extract components from solution vector
    
    U1 = U(1:nCoordinates);
    U2 = U(nCoordinates+(1:nCoordinates));
    P = U(2*nCoordinates+(1:nElements));
    
    % Generate figures
    
    plot_p_P0(P,Mesh);
    plot_u_P1(U1,U2,Mesh);
  
  elseif(strcmp(upper(type),'TH'))
      
    % Extract components from solution vector
    
    U1 = U(1:(nCoordinates+nEdges));
    U2 = U(nCoordinates+nEdges+(1:(nCoordinates+nEdges)));
    P = U(2*(nCoordinates+nEdges)+(1:nCoordinates));
    
    % Generate figures
    
    plot_p_P1(P,Mesh);
    plot_u_P2(U1,U2,Mesh);
    
  else
    
    % Extract components from solution vector  
      
    U1 = U(1:nEdges);
    U2 = U(nEdges+(1:nEdges));
    P = U(2*nEdges+(1:nElements));
    
    % Generate figures
    
    plot_p_P0(P,Mesh);  
    plot_u_CR(U1,U2,Mesh);
      
  end
  hold off;
  
  % Assign output arguments
  
  if(nargout > 0)
    varargout{1} = fig;
  end
  
return  
  
%%% Plot routine for the the MINI element %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = plot_u_MINI(U1,U2,Mesh)
% PLOT_U_MINI Plot routine for MINI elements.
%
%   FIG = PLOT_U_MINI(U1,U2,MESH) generates a plot of the velocity field
%   [U1,U2] for the MINI elements on the struct MESH and returns the handle
%   FIG to the figure.
%
%   The struct should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh, where
%                M is equal to the number of vertices contained in the
%                mesh.
%    ELEMENTS    M-by-3 matrix specifying the elements of the mesh, where M
%                is equal to the number of elements contained in the mesh.
%
%   Example:
%
%   fig = plot_u_MINI(U1,U2,Mesh);

%   Copyright 2005-2005 Patrick Meury & Kah-Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  OFFSET = 0.05;
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  
  % Compute axes limits
  
  XMax = max(Mesh.Coordinates(:,1));
  XMin = min(Mesh.Coordinates(:,1));
  YMax = max(Mesh.Coordinates(:,2));
  YMin = min(Mesh.Coordinates(:,2));
  
  XLim = [XMin XMax] + (XMax-XMin)*OFFSET*[-1 1];
  YLim = [YMin YMax] + (YMax-YMin)*OFFSET*[-1 1];
  
  % Compute barycenters of all elements
  
  XBar = 1/3*(Mesh.Coordinates(Mesh.Elements(:,1),:) + ...
              Mesh.Coordinates(Mesh.Elements(:,2),:) + ...
              Mesh.Coordinates(Mesh.Elements(:,3),:));

  % Evaluate shape functions at the barycenter of the reference element
  
  N = shap_MINI([1/3 1/3]);
          
  % Compute velocity field at barycenters
  
  u = N(1)*U1(Mesh.Elements(:,1)) + ...
      N(2)*U1(Mesh.Elements(:,2)) + ...
      N(3)*U1(Mesh.Elements(:,3)) + ...
      N(4)*U1(nCoordinates+(1:nElements));
  v = N(1)*U2(Mesh.Elements(:,1)) + ...
      N(2)*U2(Mesh.Elements(:,2)) + ...
      N(3)*U2(Mesh.Elements(:,3)) +...
      N(4)*U2(nCoordinates+(1:nElements));     
          
  % Generate figure
   
  quiver(XBar(:,1),XBar(:,2),u,v,'k-');
  set(gca,'XLim',XLim,'YLim',YLim,'DataAspectRatio',[1 1 1]);

return

%%% Plot routine for P1 elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = plot_u_P1(U1,U2,Mesh)
% PLOT_U_P1 Plot routine for P2 elements.
%
%   FIG = PLOT_U_P1(U1,U2,MESH) generates a plot of the velocity field
%   [U1,U2] for the P1 elements on the struct MESH and returns the handle
%   FIG to the figure.
%
%   The struct should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh, where
%                M is equal to the number of vertices contained in the
%                mesh.
%    ELEMENTS    M-by-3 matrix specifying the elements of the mesh, where
%                M is equal to the number of elements contained in the
%                mesh.
%    EDGES       P-by-2 matrix specifying all edges of the mesh.
%
%   Example:
%
%   fig = plot_u_P1(U1,U2,Mesh);

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  OFFSET = 0.05;
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  
  % Compute axes limits
  
  XMax = max(Mesh.Coordinates(:,1));
  XMin = min(Mesh.Coordinates(:,1));
  YMax = max(Mesh.Coordinates(:,2));
  YMin = min(Mesh.Coordinates(:,2));
  
  XLim = [XMin XMax] + (XMax-XMin)*OFFSET*[-1 1];
  YLim = [YMin YMax] + (YMax-YMin)*OFFSET*[-1 1];
  
  % Compute barycenters of all elements
  
  XBar = 1/3*(Mesh.Coordinates(Mesh.Elements(:,1),:) + ...
              Mesh.Coordinates(Mesh.Elements(:,2),:) + ...
              Mesh.Coordinates(Mesh.Elements(:,3),:));
          
  % Evaluate shape functions at the barycenter of the reference element
  
  N = shap_LFE([1/3 1/3]);
  
  % Compute velocity field at barycenters
  
  u = zeros(nElements,1);
  v = zeros(nElements,1);
  for i = 1:nElements
    idx = Mesh.Elements(i,:);
    u(i) = N(1)*U1(idx(1)) + ...
           N(2)*U1(idx(2)) + ...
           N(3)*U1(idx(3));
    v(i) = N(1)*U2(idx(1)) + ...
           N(2)*U2(idx(2)) + ...
           N(3)*U2(idx(3));
  end
          
  % Generate figure
   
  quiver(XBar(:,1),XBar(:,2),u,v,'k-');
  set(gca,'XLim',XLim,'YLim',YLim,'DataAspectRatio',[1 1 1]);

return

%%% Plot routine for P2 elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = plot_u_P2(U1,U2,Mesh)
% PLOT_U_P2 Plot routine for P2 elements.
%
%   FIG = PLOT_U_P2(U1,U2,MESH) generates a plot of the velocity field
%   [U1,U2] for the P2 elements on the struct MESH and returns the handle
%   FIG to the figure.
%
%   The struct should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh, where
%                M is equal to the number of vertices contained in the
%                mesh.
%    ELEMENTS    M-by-3 matrix specifying the elements of the mesh, where
%                M is equal to the number of elements contained in the
%                mesh.
%    EDGES       P-by-2 matrix specifying all edges of the mesh.
%    VERT2EDGE   M-by-M sparse matrix which specifies wheter the two
%                vertices i and j are connected by an edge with number
%                VERT2EDGE(i,j).
%
%   Example:
%
%   fig = plot_u_P2(U1,U2,Mesh);

%   Copyright 2005-2005 Patrick Meury & Kah-Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  OFFSET = 0.05;
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  nEdges = size(Mesh.Edges,1);
  
  % Compute axes limits
  
  XMax = max(Mesh.Coordinates(:,1));
  XMin = min(Mesh.Coordinates(:,1));
  YMax = max(Mesh.Coordinates(:,2));
  YMin = min(Mesh.Coordinates(:,2));
  
  XLim = [XMin XMax] + (XMax-XMin)*OFFSET*[-1 1];
  YLim = [YMin YMax] + (YMax-YMin)*OFFSET*[-1 1];
  
  % Compute barycenters of all elements
  
  XBar = 1/3*(Mesh.Coordinates(Mesh.Elements(:,1),:) + ...
              Mesh.Coordinates(Mesh.Elements(:,2),:) + ...
              Mesh.Coordinates(Mesh.Elements(:,3),:));
          
  % Evaluate shape functions at the barycenter of the reference element
  
  N = shap_QFE([1/3 1/3]);
  
  % Compute velocity field at barycenters
  
  u = zeros(nElements,1);
  v = zeros(nElements,1);
  for i = 1:nElements
    vidx = Mesh.Elements(i,:);
    idx = [vidx ...
           Mesh.Vert2Edge(vidx(1),vidx(2))+nCoordinates ...
           Mesh.Vert2Edge(vidx(2),vidx(3))+nCoordinates ...
           Mesh.Vert2Edge(vidx(3),vidx(1))+nCoordinates];
    u(i) = N(1)*U1(idx(1)) + ...
           N(2)*U1(idx(2)) + ...
           N(3)*U1(idx(3)) + ... 
           N(4)*U1(idx(4)) + ...
           N(5)*U1(idx(5)) + ...
           N(6)*U1(idx(6));
    v(i) = N(1)*U2(idx(1)) + ...
           N(2)*U2(idx(2)) + ...
           N(3)*U2(idx(3)) + ... 
           N(4)*U2(idx(4)) + ...
           N(5)*U2(idx(5)) + ...
           N(6)*U2(idx(6));
  end
       
  % Generate figure
   
  quiver(XBar(:,1),XBar(:,2),u,v,'b-');
  set(gca,'XLim',XLim,'YLim',YLim,'DataAspectRatio',[1 1 1]);

return

%%% Plot routine for Crouzeix-Raviart elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = plot_u_CR(U1,U2,Mesh)
% PLOT_U_CR Plot routine for Crouzeix-Raviart elements.
%
%   FIG = PLOT_U_CR(U1,U2,MESH) generates a plot of the velocity field
%   [U1,U2] for the Crouzeix-Raviart elements on the struct MESH and
%   returns the handle FIG to the figure.
%
%   The struct should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh, where
%                M is equal to the number of vertices contained in the
%                mesh.
%    ELEMENTS    M-by-3 matrix specifying the elements of the mesh, where
%                M is equal to the number of elements contained in the
%                mesh.
%    EDGES       P-by-2 matrix specifying all edges of the mesh.
%    VERT2EDGE   M-by-M sparse matrix which specifies wheter the two
%                vertices i and j are connected by an edge with number
%                VERT2EDGE(i,j).
%
%   Example:
%
%   fig = plot_u_CR(U1,U2,Mesh);

%   Copyright 2005-2005 Patrick Meury & Kah-Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  OFFSET = 0.05;  
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  nEdges = size(Mesh.Edges,1);
  
  % Compute axes limits
  
  XMax = max(Mesh.Coordinates(:,1));
  XMin = min(Mesh.Coordinates(:,1));
  YMax = max(Mesh.Coordinates(:,2));
  YMin = min(Mesh.Coordinates(:,2));
  
  XLim = [XMin XMax] + (XMax-XMin)*OFFSET*[-1 1];
  YLim = [YMin YMax] + (YMax-YMin)*OFFSET*[-1 1];
  
  % Compute barycenters of all elements
  
  XBar = 1/3*(Mesh.Coordinates(Mesh.Elements(:,1),:) + ...
              Mesh.Coordinates(Mesh.Elements(:,2),:) + ...
              Mesh.Coordinates(Mesh.Elements(:,3),:));

  % Evaluate shape functions at the barycenter of the reference element
  
  N = shap_CR([1/3 1/3]);
  
  % Compute velocity field at barycenters
  
  u = zeros(nElements,1);
  v = zeros(nElements,1);
  for i = 1:nElements
    vidx = Mesh.Elements(i,:);
    idx = [Mesh.Vert2Edge(vidx(2),vidx(3)) ...
           Mesh.Vert2Edge(vidx(3),vidx(1)) ...
           Mesh.Vert2Edge(vidx(1),vidx(2))];
    u(i) = N(1)*U1(idx(1)) + ...
           N(2)*U1(idx(2)) + ...
           N(3)*U1(idx(3));
    v(i) = N(1)*U2(idx(1)) + ...
           N(2)*U2(idx(2)) + ...
           N(3)*U2(idx(3));
  end
  
  % Generate figure
   
  quiver(XBar(:,1),XBar(:,2),u,v,'k-');
  set(gca,'XLim',XLim,'YLim',YLim,'DataAspectRatio',[1 1 1]);
          
return

%%% Plot routine for P0 elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = plot_p_P0(P,Mesh)
% PLOT_P_P0 Plot routine for P0 elements.
%
%   FIG = PLOT_P_P0(P,MESH) generates a plot of the pressure P for the P0
%   elements on the struct MESH and returns the handle FIG to the figure.
%
%   The struct should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh, where
%                M is equal to the number of vertices contained in the
%                mesh.
%    ELEMENTS    M-by-3 matrix specifying the elements of the mesh, where
%                M is equal to the number of elements contained in the
%                mesh.
%
%   Example:
%
%   fig = plot_p_P0(P,Mesh);

%   Copyright 2005-2005 Patrick Meury & Kah-Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  OFFSET = 0.05;

  % Compute axes limits
  
  XMax = max(Mesh.Coordinates(:,1));
  XMin = min(Mesh.Coordinates(:,1));
  YMax = max(Mesh.Coordinates(:,2));
  YMin = min(Mesh.Coordinates(:,2));
  PMax = max(P);
  PMin = min(P);
  XLim = [XMin XMax] + (XMax-XMin)*OFFSET*[-1 1];
  YLim = [YMin YMax] + (YMax-YMin)*OFFSET*[-1 1];
  if(PMin < PMax)
    CLim = [PMin PMax] + (PMax-PMin)*OFFSET*[-1 1];
  else
    CLim = PMin*[1-OFFSET 1+OFFSET];  
  end
  
  % Generate figure

  patch('Vertices',Mesh.Coordinates, ...
        'Faces',Mesh.Elements, ...
        'CData',P, ...
        'EdgeColor','none', ...
        'FaceColor','flat');
  set(gca,'XLim',XLim,'YLim',YLim,'CLim',CLim,'DataAspectRatio',[1 1 1]);
  
return

%%% Plot routine for P1 elements %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = plot_p_P1(P,Mesh)
% PLOT_P_P1 Plot routine for P1 elements.
%
%   FIG = PLOT_P_P1(P,MESH) generates a plot of the pressure P for the P1
%   elements on the struct MESH and returns the handle FIG to the figure.
%
%   The struct should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh, where
%                M is equal to the number of vertices contained in the
%                mesh.
%    ELEMENTS    M-by-3 matrix specifying the elements of the mesh, where
%                M is equal to the number of elements contained in the
%                mesh.
%
%   Example:
%
%   fig = plot_p_P1(P,Mesh);

%   Copyright 2005-2005 Patrick Meury & Kah-Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  OFFSET = 0.05;

  % Compute axes limits
  
  XMax = max(Mesh.Coordinates(:,1));
  XMin = min(Mesh.Coordinates(:,1));
  YMax = max(Mesh.Coordinates(:,2));
  YMin = min(Mesh.Coordinates(:,2));
  PMax = max(P);
  PMin = min(P);
  XLim = [XMin XMax] + (XMax-XMin)*OFFSET*[-1 1];
  YLim = [YMin YMax] + (YMax-YMin)*OFFSET*[-1 1];
  if(PMin < PMax)
    CLim = [PMin PMax] + (PMax-PMin)*OFFSET*[-1 1];
  else
    CLim = PMin*[1-OFFSET 1+OFFSET];  
  end
  
  % Generate figure

  patch('Vertices',Mesh.Coordinates, ...
        'Faces',Mesh.Elements, ...
        'CData',P, ...
        'EdgeColor','none', ...
        'FaceColor','interp');
  set(gca,'XLim',XLim,'YLim',YLim,'CLim',CLim,'DataAspectRatio',[1 1 1]);
  
return
