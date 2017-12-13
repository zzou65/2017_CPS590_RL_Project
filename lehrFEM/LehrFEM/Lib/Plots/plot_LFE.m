function varargout = plot_LFE(U,Mesh,fig)
% PLOT_LFE Plot finite element solution.
%
%   PLOT_LFE(U,MESH) generates a plot of the finite element solution U on
%   the mesh MESH.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 matrix specifying the elements of the mesh.
%
%   H = PLOT_LFE(U,MESH) also returns the handle to the figure.
%
%   Example:
%
%   plot_LFE(U,MESH);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  OFFSET = 0.00;
  
  % Compute axes limits
  
  XMin = min(Mesh.Coordinates(:,1));
  XMax = max(Mesh.Coordinates(:,1));
  YMin = min(Mesh.Coordinates(:,2));
  YMax = max(Mesh.Coordinates(:,2));
  XLim = [XMin XMax] + OFFSET*(XMax-XMin)*[-1 1]; 
  YLim = [YMin YMax] + OFFSET*(YMax-YMin)*[-1 1];
  
  % Generate figure
    
  if(isreal(U))
  
    % Compute color axes limits 
      
    CMin = min(U);
    CMax = max(U);
    if(CMin < CMax)          % or error will occur in set function
      CLim = [CMin CMax] + OFFSET*(CMax-CMin)*[-1 1];
    else
      CLim = [1-OFFSET 1+OFFSET]*CMin;   
    end
         
    % Plot real finite element solution  
    % Create new figure, if argument 'fig' is not specifiied
    % Otherwise this argument is supposed to be a figure handle
    if (nargin < 3),  fig = figure('Name','Linear finite elements');
    else figure(fig); end
    
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
    CLim = [CMin CMax] + OFFSET*(CMax-CMin)*[-1 1];
    
    % Plot imaginary finite element solution  
      
    fig_1 = figure('Name','Linear finite elements');
    patch('faces', Mesh.Elements, ...
          'vertices', [Mesh.Coordinates(:,1) Mesh.Coordinates(:,2) real(U)], ...
          'CData', real(U), ...
          'facecolor', 'interp', ...
          'edgecolor', 'none');
    set(gca,'XLim',XLim,'YLim',YLim,'CLim',CLim,'DataAspectRatio',[1 1 1]);
    fig_2 = figure('Name','Linear finite elements');
    patch('faces', Mesh.Elements, ...
          'vertices', [Mesh.Coordinates(:,1) Mesh.Coordinates(:,2) imag(U)], ...
          'CData', imag(U), ...
          'facecolor', 'interp', ...
          'edgecolor', 'none');  
    set(gca,'XLim',XLim,'YLim',YLim,'CLim',CLim,'DataAspectRatio',[1 1 1]);
    if(nargout > 0)
      varargout{1} = fig_1;
      varargout{2} = fig_2;
    end
      
  end
  
return