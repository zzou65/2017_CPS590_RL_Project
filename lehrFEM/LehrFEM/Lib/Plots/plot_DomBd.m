function varargout = plot_DomBd(Mesh,varargin)
% PLOT_DomBd mesh domain boundary plot.
%
%   PLOT_DOMBD(MESH) generate 2D plot of the mesh boundary.
%
%   PLOT_DOMBD(MESH,OPT) adds labels to the plot, where OPT is a character 
%   string made from one element from any or all of the following 
%   characters:
%    a Dipslay axes on the plot.
%    s Add title and axes labels to the plot.
%
%   H = PLOT_DOMBD(MESH,OPT) also returns the handle to the figure.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 or N-by-4 matrix specifying the elements of the
%                mesh. 
%
%   Example:
%
%   plot_DomBd(Mesh,'as');
%
%   See also get_BdEdges, add_Edges.

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize constants
   
  OFFSET = 0.05;      % Offset parameter 
  BDEDGECOLOR = 'r';  % Boundary edge color
  
  % Check mesh data structure and add necessary fields
  
  if(~isfield(Mesh,'Edges'))
    Mesh = add_Edges(Mesh);
  end
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  nEdges = size(Mesh.Edges,1);
  
  % Compute axes limits
  
  X = Mesh.Coordinates(:,1);
  Y = Mesh.Coordinates(:,2);
  XMin = min(X);
  XMax = max(X);
  YMin = min(Y);
  YMax = max(Y);
  XLim = [XMin XMax] + OFFSET*(XMax-XMin)*[-1 1];
  YLim = [YMin YMax] + OFFSET*(YMax-YMin)*[-1 1];
        
  % Compute boundary edges for piecewise linear boundaries
     
  Loc = get_BdEdges(Mesh);
  BdEdges_x = zeros(2,size(Loc,1));
  BdEdges_y = zeros(2,size(Loc,1));
  BdEdges_x(1,:) = Mesh.Coordinates(Mesh.Edges(Loc,1),1)';
  BdEdges_x(2,:) = Mesh.Coordinates(Mesh.Edges(Loc,2),1)';
  BdEdges_y(1,:) = Mesh.Coordinates(Mesh.Edges(Loc,1),2)';
  BdEdges_y(2,:) = Mesh.Coordinates(Mesh.Edges(Loc,2),2)';

  % Generate plot

  if(~ishold)
      hold on;
  end
  fig = plot(BdEdges_x,BdEdges_y,[BDEDGECOLOR '-']);
  hold off;
  set(gca,'XLim',XLim, ...
      'YLim',YLim, ...
      'DataAspectRatio',[1 1 1], ...
      'Box','on', ...
      'Visible','off');
  
  if(nargin > 1)
      opt = varargin{1};
      if(~isempty(findstr('a',opt)))
          set(gca,'Visible','on');
          if(~isempty(findstr('s',opt)))
              if(size(Mesh.Elements,2) == 3)
                  title(['{\bf 2D triangular mesh}']);
              else
                  title(['{\bf 2D quadrilateral mesh}']);
              end
              xlabel(['{\bf # Vertices  :  ', int2str(nCoordinates), ...
                  ',      # Elements  :  ', int2str(nElements), ...
                  ',      # Edges  :  ',int2str(nEdges),'}']);
          end
      end
  end
  drawnow;

  % Assign output arguments

  if(nargout > 0)
      varargout{1} = fig;
  end

return