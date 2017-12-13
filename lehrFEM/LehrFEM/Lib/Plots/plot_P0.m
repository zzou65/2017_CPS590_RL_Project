function varargout = plot_P0(U,Mesh)
% PLOT_P0 Plot finite element solution.
%
%   PLOT_P0(U,MESH) generates a plot of the finite element solution U on
%   the mesh MESH.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 or N-by-4 matrix specifying the elements of the
%                mesh.
%
%   H = PLOT_P0(U,MESH) also returns the handle to the figure.
%
%   Example:
%
%   plot_P0(U,MESH);

%   Copyright 2006-2006 Patrick Meury & Kah-Ling Sia
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
  UMax = max(U);
  UMin = min(U);
  XLim = [XMin XMax] + (XMax-XMin)*OFFSET*[-1 1];
  YLim = [YMin YMax] + (YMax-YMin)*OFFSET*[-1 1];
  if(UMin < UMax)
    CLim = [UMin UMax] + (UMax-UMin)*OFFSET*[-1 1];
  else
    CLim = UMin*[1-OFFSET 1+OFFSET];  
  end
  
  % Generate figure

  fig = figure('Name','Piecewise constant finite elements');
  patch('Vertices',Mesh.Coordinates, ...
        'Faces',Mesh.Elements, ...
        'CData',U, ...
        'EdgeColor','none', ...
        'FaceColor','flat');
  set(gca,'XLim',XLim,'YLim',YLim,'CLim',CLim,'DataAspectRatio',[1 1 1]);

  if(nargout > 0)
    varargout{1} = fig;  
  end
  
return
