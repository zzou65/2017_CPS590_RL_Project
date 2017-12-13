function varargout = plot_PBD(U,Mesh)
% PLOT_PBD Plot finite element solution.
%
%   PLOT_PBD(U,MESH) generates a plot fo the finite element solution U on
%   the mesh MESH.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%    EDGES        P-by-2 matrix specifying all edges of the mesh.
%    VERT2EDGE    M-by-M sparse matrix which specifies wheter the two
%                 vertices i and j are connected by an edge with number
%                 VERT2EDGE(i,j).
%    EDGE2ELEM    N-by-2 matrix connecting edges to elements. The first
%                 column specifies the left hand side element where the
%                 second column specifies the right hand side element.
%    EDGELOC      P-by-2 matrix connecting egdes to local edges of
%                 elements.
%    DELTA        P-by-1 matrix specifying the boundary correction term on
%                 every edge.
%
%   H = PLOT_PBD(U,MESH) also returns the handle to the figure.
%
%   Example:
%
%   plot_PBD(U,Mesh);
%
%   See also refine_REG, add_Edge2Elem, get_BdEdges.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

 % Initialize constants
  
  Rot = [0 -1; 1 0];
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  Offset = 0.05;
  
  % Refine the mesh for plotting
  
  Delta = Mesh.Delta; 
  Mesh = add_Edge2Elem(Mesh);
  BdEdges = get_BdEdges(Mesh);
  [Rows,Cols,Vals] = find(Mesh.Edge2Elem(BdEdges,:));
  Vals = sortrows([Rows Cols Vals],1);
  BdElem = Vals(:,3);
  normal = Mesh.Coordinates(Mesh.Elements(BdElem,3),:) - ...
           Mesh.Coordinates(Mesh.Elements(BdElem,2),:);         
  normal = (normal*Rot)./(sqrt(sum(normal.^2,2))*ones(1,2));
  Mesh = refine_REG(Mesh);
  MidPoints = BdEdges + nCoordinates;
  Mesh.Coordinates(MidPoints,:) = Mesh.Coordinates(MidPoints,:) + ...
                                  Delta(BdEdges)*ones(1,2).*normal;
  
  % Compute axes limits
  
  XMin = min(Mesh.Coordinates(:,1));
  XMax = max(Mesh.Coordinates(:,1));
  YMin = min(Mesh.Coordinates(:,2));
  YMax = max(Mesh.Coordinates(:,2));
  XLim = [XMin XMax] + Offset*(XMax-XMin)*[-1 1]; 
  YLim = [YMin YMax] + Offset*(YMax-YMin)*[-1 1];
  
  % Generate figure
    
  if(isreal(U))
  
    % Compute color axes limits 
      
    CMin = min(U);
    CMax = max(U);
    if(CMin < CMax)
      CLim = [CMin CMax] + Offset*(CMax-CMin)*[-1 1];
    else
      CLim = [1-Offset 1+Offset]*CMin;   
    end
    
    % Plot real finite element solution  
      
    h = figure('Name', ...
               'Quadratic finite elements with parabolic boundaries');
    patch('faces', Mesh.Elements, ...
          'vertices', [Mesh.Coordinates(:,1) Mesh.Coordinates(:,2) U], ...
          'CData', U, ...
          'facecolor', 'interp', ...
          'edgecolor', 'none');
    set(gca,'XLim',XLim,'YLim',YLim,'CLim',CLim,'DataAspectRatio',[1 1 1]);
    
    if(nargout > 0)
      varargout{1} = fig;
    end
 
  else  
      
    % Compute color axes limits 
      
    CMin = min([real(U); imag(U)]);
    CMax = max([real(U); imag(U)]);
    CLim = [CMin CMax] + Offset*(CMax-CMin)*[-1 1];
    
    % Plot imaginary finite element solution  
      
    h1 = figure('Name', ...
                'Quadratic finit elements with parabolic boundaries');
    patch('faces', Mesh.Elements, ...
          'vertices', [Mesh.Coordinates(:,1) Mesh.Coordinates(:,2) real(U)], ...
          'CData', real(U), ...
          'facecolor', 'interp', ...
          'edgecolor', 'none');
    set(gca,'XLim',XLim,'YLim',YLim,'CLim',CLim,'DataAspectRatio',[1 1 1]);
    h2 = figure('Name', ...
                'Quadratic finit elements with parabolic boundaries');
    patch('faces', Mesh.Elements, ...
          'vertices', [Mesh.Coordinates(:,1) Mesh.Coordinates(:,2) imag(U)], ...
          'CData', imag(U), ...
          'facecolor', 'interp', ...
          'edgecolor', 'none');  
    set(gca,'XLim',XLim,'YLim',YLim,'CLim',CLim,'DataAspectRatio',[1 1 1]);
    if(nargout > 0)
      varargout{1} = h1;
      varargout{2} = h2;
    end
      
  end
  
return