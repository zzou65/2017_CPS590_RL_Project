function varargout = plot_EST(Eta,Mesh)
% PLOT_EST Plot routine for error estimators.
%
%   FIG = PLOT_EST(ETA,MESH) generates a plot of the error estimators ETA 
%   on the struct MESH.
%
%   The struct should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    M-by-3 matrix specifying the elements of the mesh.
%
%   FIG = PLOT_EST(ETA,MESH) also returns the handle FIG to the figure.
%
%   Example:
%
%   fig = plot_EST(Eta,Mesh);

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
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
  EtaMax = max(Eta);
  EtaMin = min(Eta);
  XLim = [XMin XMax] + (XMax-XMin)*OFFSET*[-1 1];
  YLim = [YMin YMax] + (YMax-YMin)*OFFSET*[-1 1];
  if(EtaMin < EtaMax)
    CLim = [EtaMin EtaMax] + (EtaMax-EtaMin)*OFFSET*[-1 1];
  else
    CLim = EtaMin*[1-OFFSET 1+OFFSET];  
  end
  
  % Generate figure

  fig = figure;
  patch('Vertices',Mesh.Coordinates, ...
        'Faces',Mesh.Elements, ...
        'CData',Eta, ...
        'EdgeColor','none', ...
        'FaceColor','flat');
  set(gca,'XLim',XLim,'YLim',YLim,'CLim',CLim,'DataAspectRatio',[1 1 1]);
  
  % Assign output arguments
  
  if(nargin > 1)
    varargout{1} = fig;  
  end
  
return