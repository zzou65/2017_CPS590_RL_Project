function varargout = plot_Qual(Mesh)
% PLOT_QUAL Element quality plot.
%
%   PLOT_QUAL(MESH) generates 2D plot of the element quality of all elements
%   in the struct MESH. The triangle quality of the mesh is computed according
%   to the formula
%
%     q(T) = 2*R_in/R_out,
%
%   where R_in denotes the radius of the largest inscribed circle and R_out
%   the radius of the smallest circumscribed circle of the triangle T.
%
%   H = PLOT_QUAL(MESH) also returns the handle to the figure.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh, where M
%                is equal to the number of vertices contained in the mesh.
%    ELEMENTS    M-by-3 matrix specifying the elements of the mesh, where M
%                is equal to the number of elements contained in the mesh. 
%
%   Example:
%
%   plot_Qual(Mesh);
%
%   See also TRIQUAL, QUADQUAL.

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constant
   
  OFFSET = 0.05;  % Offset parameter 
    
  % Compute axes limits
  
  X = Mesh.Coordinates(:,1);
  Y = Mesh.Coordinates(:,2);
  XMin = min(X);
  XMax = max(X);
  YMin = min(Y);
  YMax = max(Y);
  XLim = [XMin XMax] + OFFSET*(XMax-XMin)*[-1 1];
  YLim = [YMin YMax] + OFFSET*(YMax-YMin)*[-1 1];
      
  % Compute element qualities
  
  [nElements,nVert] = size(Mesh.Elements);
  Q = zeros(nElements,1);
  if(nVert == 3)
    for i = 1:nElements
      Q(i) = TriQual(Mesh.Coordinates(Mesh.Elements(i,:),:));    
    end
  else
    for i = 1:nElements
      Q(i) = QuadQual(Mesh.Coordinates(Mesh.Elements(i,:),:));    
    end
  end  
      
  % Plot element qualities
  
  fig = figure('Name','Element quality');
  patch('faces',Mesh.Elements, ...
        'vertices',Mesh.Coordinates, ...
        'facevertexcdata',Q, ...
        'facecolor','flat', ...
        'edgecolor','black');
  axis('equal');
  box('on');
  set(gca,'XLim',XLim,'YLim',YLim,'ZLim',[0 1],'CLim',[0 1]);
  colorbar;
  title('{\bf Element quality}');
  xlabel(['{\bf Mean  :  ',num2str(mean(Q),'%1.2f'), ...
          ',      Min  :  ',num2str(min(Q),'%1.2f'), ...
          ',      Max  :  ',num2str(max(Q),'%1.2f'),'}']);
  drawnow;
   
  % Assign output arguments    
          
  if(nargout > 0)
    varargout{1} = fig;
  end
  
return  
