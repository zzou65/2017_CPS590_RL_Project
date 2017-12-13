function varargout = plot_Meshes(Meshcell,varargin)
% plot_Meshes plots several Meshes in one figure  
% 
%   PLOT_MESHES(MESHCELL) generate 2D plot of the meshes.
%   
%
%   PLOT_MESHES(MESHCELL,OPT) adds labels to the plot, where OPT is a character string
%   made from one element from any or all of the following characters:
%    p Add vertex labels to the plot.
%    e Add edge labels/flags to the plot.
%    t Add element labels/flags to the plot.
%    a Dipslay axes on the plot.
%    s Add title and axes labels to the plot.
%    f Do NOT create new window for the mesh plot
%    [c add patch color to elements according to their flags] TODO !
%
%   H = PLOT_MESHES(MESHCELL,OPT) also returns the handle to the figure.
%
%   Meshcell is a cell of meshes
%   The structs MESHCELL{i} should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 or N-by-4 matrix specifying the elements of the
%                mesh. 
%
%   Example:
%
%   plot_Meshes(Mesh,'petas');
%
%   See also get_BdEdges, add_Edges.%
%  Meshcell is a cell of  Meshes
%
%   Copyright 2005-2009 Patrick Meury, Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

nMeshes=max(size(Meshcell));

Mesh=Meshcell{1};

if(nargin > 1)
  opt = varargin{1};
else
  opt = ' ';
end

  % Initialize constants
   
  OFFSET = 0.05;      % Offset parameter 
  EDGECOLOR = 'b';    % Interior edge color
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
 
  if(isempty(findstr('f',opt)))
    fig = figure('Name','Mesh plot'); 
  end
  
  if(~ishold)
    hold on;  
  end
  patch('Faces', Mesh.Elements, ...
        'Vertices', Mesh.Coordinates, ...
        'FaceColor', 'none', ...
        'EdgeColor', EDGECOLOR);
    for num=2:nMeshes
        patch('Faces', Meshcell{num}.Elements, ...
            'Vertices', Meshcell{num}.Coordinates, ...
            'FaceColor', 'none', ...
            'EdgeColor', 'g');
    end
  plot(BdEdges_x,BdEdges_y,[BDEDGECOLOR '-']);
  if nargin>2
      trans_bary=varargin{2};
      plot(trans_bary(:,1),trans_bary(:,2),'r*');
  end
  hold off;
  set(gca,'XLim',XLim, ...
          'YLim',YLim, ...
          'DataAspectRatio',[1 1 1], ...
          'Box','on', ...
          'Visible','off');
  
  
  % Add labels/flags according to the string OPT
  
    
    % Add vertex labels
    
    if(~isempty(findstr('p',opt)))
      add_VertLabels(Mesh.Coordinates);  
    end
    
    % Add element labels/flags to the plot
    
    if(~isempty(findstr('t',opt)))
      %if(isfield(Mesh,'ElemFlag'))
       % add_ElemLabels(Mesh.Coordinates,Mesh.Elements,Mesh.ElemFlag);  
      %else
        add_ElemLabels(Mesh.Coordinates,Mesh.Elements,1:nElements); 
      %end
    end
    
    % Add edge labels/flags to the plot
    
    if(~isempty(findstr('e',opt)))
      if(isfield(Mesh,'BdFlags'))
        %add_EdgeLabels(Mesh.Coordinates,Mesh.Edges,Mesh.BdFlags);
        add_EdgeLabels(Mesh.Coordinates,Mesh.Edges,1:nEdges);
      else
        add_EdgeLabels(Mesh.Coordinates,Mesh.Edges,1:nEdges);    
      end
    end
    
    % Turn on axes, titles and labels
    
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
    
  drawnow;
  
  % Assign output arguments    
          
  if(nargout > 0)
    varargout{1} = fig;
  end
  
return  


%%% Add vertex labels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = add_VertLabels(Coordinates)
% ADD_VERTLABELS Add vertex labels to the plot.
%
%   ADD_VERTLABELS(COORDINATES) adds vertex labels to the current
%   figure.
%
%   Example:
%
%   add_VertLabels(Mesh.Coordinates);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants

  WEIGHT = 'bold';
  SIZE = 8;
  COLOR = 'k';
  
  % Add vertex labels to the plot
  
  nCoordinates = size(Coordinates,1);
  for i = 1:nCoordinates
    text(Coordinates(i,1),Coordinates(i,2),int2str(i), ...
         'HorizontalAlignment','Center', ...
         'VerticalAlignment','Middle', ...
         'Color',COLOR, ...
         'FontWeight',WEIGHT, ...
         'FontSize',SIZE);
  end    
      
return

%%% Add element labels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = add_ElemLabels(Coordinates,Elements,Labels)
% ADD_ELEMLABELS Add element labels to the plot.
%
%   ADD_ELEMLABELS(COORDINATES,ELEMENTS,LABELS) adds the element labels
%   LABELS to the current figure.
%
%   Example:
%
%   add_ElemLabels(Mesh.Coordinates,Mesh.Elements,Labels);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize constants

  WEIGHT = 'bold'; 
  SIZE = 8;        
  COLOR = 'k';

  % Add element labels to the plot
  
  [nElements,nVert] = size(Elements);
  for i = 1:nElements
    CoordMid = sum(Coordinates(Elements(i,:),:),1)/nVert;
    text(CoordMid(1),CoordMid(2),int2str(Labels(i)), ...
         'HorizontalAlignment','Center', ...
         'VerticalAlignment','Middle', ...
         'Color',COLOR, ...
         'FontWeight',WEIGHT, ...
         'FontSize',SIZE);
  end
  
return

%%% Add edge labels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [] = add_EdgeLabels(Coordinates,Edges,Labels)
% ADD_EDGELABELS Add edge labels to the plot.
%
%   ADD_EDGELABELS(COORDINATES,EDGES,LABELS) adds the edge labels LABELS to
%   the current figure. 
%
%   Example:
%
%   add_EdgeLabels(Coordinates,Edges,Labels);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  WEIGHT = 'bold';
  SIZE = 8;
  COLOR = 'k';

  % Add edge labels to the plot
 
  nEdges = size(Edges,1);
  for i = 1:nEdges
    CoordMid = (Coordinates(Edges(i,1),:)+Coordinates(Edges(i,2),:))/2;  
    text(CoordMid(1),CoordMid(2),int2str(Labels(i)), ...
         'HorizontalAlignment','Center', ...
         'VerticalAlignment','Middle', ...
         'Color',COLOR, ...
         'FontWeight',WEIGHT, ...
         'FontSize',SIZE);
  end
    
return