function vargout = plot_MeshPatch(Mesh,pbVertices,plot_param,varargin)

% plots the mesh and the deformed mesh and the barycenters of the patches

%   Copyright 2008 Holger Heumann, Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

tol = 10^-4;

Mesh = init_LEB(Mesh);
Mesh = add_Patches(Mesh);
Mesh = add_Edge2Elem(Mesh);
pMesh = Mesh;
pMesh.Coordinates = pbVertices(:,[1 2]);


intersec = aff_elems2(Mesh, pMesh,pbVertices);


% Initialize constants
   
OFFSET = 0.05;      % Offset parameter 
EDGECOLOR = 'b';    % Interior edge color
BDEDGECOLOR = 'r';  % Boundary edge color
  
% Check mesh data structure and add necessary fields
  
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
 
if(isempty(findstr('f',plot_param)))
  fig = figure('Name','Mesh plot'); 
end
  
if(~ishold)
  hold on;  
end
patch('Faces', Mesh.Elements, ...
      'Vertices', Mesh.Coordinates, ...
        'FaceColor', 'none', ...
        'EdgeColor', EDGECOLOR);

patch('Faces', pMesh.Elements, ...
      'Vertices', pMesh.Coordinates, ...
      'FaceColor', 'none', ...
      'EdgeColor', 'g');

plot(BdEdges_x,BdEdges_y,[BDEDGECOLOR '-']);

for i=1:nElements
    nPatches = intersec(i).nElems;
    Neigh = Mesh.Neigh(i,:);
    polygons = patches(Mesh,pMesh,i,intersec(i));
    for j=1:nPatches
        polygon = polygons(j).Polygon;
        elementid = polygons(j).Elem;
        baryc = sum(polygon,1)/size(polygon,1);
        if polyarea(polygon(:,1),polygon(:,2))>tol
            if sum(elementid==[Neigh,i])
                plot(baryc(1),baryc(2),'*');
            else
                plot(baryc(1),baryc(2),'r*');
            end;
        end
    end
end



hold off
  
set(gca,'XLim',XLim, ...
          'YLim',YLim, ...
          'DataAspectRatio',[1 1 1], ...
          'Box','on', ...
          'Visible','off');
      
% Add labels/flags according to the string plot_param
  
    
    % Add vertex labels
    
    if(~isempty(findstr('p',plot_param)))
      add_VertLabels(Mesh.Coordinates);  
    end
    
    % Add element labels/flags to the plot
    
    if(~isempty(findstr('t',plot_param)))
      %if(isfield(Mesh,'ElemFlag'))
       % add_ElemLabels(Mesh.Coordinates,Mesh.Elements,Mesh.ElemFlag);  
      %else
        add_ElemLabels(Mesh.Coordinates,Mesh.Elements,1:nElements); 
      %end
    end
    
    % Add edge labels/flags to the plot
    
    if(~isempty(findstr('e',plot_param)))
      if(isfield(Mesh,'BdFlags'))
        %add_EdgeLabels(Mesh.Coordinates,Mesh.Edges,Mesh.BdFlags);
        add_EdgeLabels(Mesh.Coordinates,Mesh.Edges,1:nEdges);
      else
        add_EdgeLabels(Mesh.Coordinates,Mesh.Edges,1:nEdges);    
      end
    end
    
    % Turn on axes, titles and labels
    
    if(~isempty(findstr('a',plot_param)))
      set(gca,'Visible','on');
      if(~isempty(findstr('s',plot_param)))
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