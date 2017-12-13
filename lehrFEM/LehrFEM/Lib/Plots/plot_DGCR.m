function varargout = plot_DGCR(U,Mesh)
% PLOT_DGCR Plot finite element solution.
%
%   PLOT_DGCR(U,MESH) generates a plot of the finite element solution U on
%   the mesh MESH.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 matrix specifying the elements of the mesh.
%
%   H = PLOT_DGCR(U,MESH) also returns the handle to the figure.
%
%   Example:
%
%   plot_DGCR(U,MESH);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  OFFSET = 0.05;
  nElements = size(Mesh.Elements,1);
  
  % Compute axes limits
  
  XMin = min(Mesh.Coordinates(:,1));
  XMax = max(Mesh.Coordinates(:,1));
  YMin = min(Mesh.Coordinates(:,2));
  YMax = max(Mesh.Coordinates(:,2));
  XLim = [XMin XMax] + OFFSET*(XMax-XMin)*[-1 1]; 
  YLim = [YMin YMax] + OFFSET*(YMax-YMin)*[-1 1];
  
  % Generate auxiliary mesh
  
  V = zeros(3*nElements,1);
  Coordinates = zeros(3*nElements,2);
  Elements = zeros(nElements,3);
  for i =1:nElements
    
     % Generate new vertex numbers 
      
     j1 = 3*(i-1)+1;
     j2 = 3*(i-1)+2;
     j3 = 3*(i-1)+3;
      
     % Update auxiliary elements
     
     Elements(i,1) = j1;
     Elements(i,2) = j2;
     Elements(i,3) = j3;
     
     % Update auxiliary vertices
     
     Coordinates(j1,:) = Mesh.Coordinates(Mesh.Elements(i,1),:);
     Coordinates(j2,:) = Mesh.Coordinates(Mesh.Elements(i,2),:);
     Coordinates(j3,:) = Mesh.Coordinates(Mesh.Elements(i,3),:);
      
     % Update the FE solution
     
     V(j1) = U(j2)+U(j3)-U(j1);
     V(j2) = U(j1)+U(j3)-U(j2);
     V(j3) = U(j1)+U(j2)-U(j3);
     
  end
  
  % Generate figure
    
  if(isreal(V))
  
    % Compute color axes limits 
      
    CMin = min(V);
    CMax = max(V);
    if(CMin < CMax)         
      CLim = [CMin CMax] + OFFSET*(CMax-CMin)*[-1 1];
    else
      CLim = [1-OFFSET 1+OFFSET]*CMin;   
    end
         
    % Plot real finite element solution  
      
    fig = figure('Name','Crouzeix-Raviart finite elements');
    patch('faces', Elements, ...
          'vertices', [Coordinates(:,1) Coordinates(:,2) V], ...
          'CData', V, ...
          'facecolor', 'interp', ...
          'edgecolor', 'none');
    set(gca,'XLim',XLim,'YLim',YLim,'CLim',CLim,'DataAspectRatio',[1 1 1]);
    
    if(nargout > 0)
      varargout{1} = fig;
    end
 
  else  
      
    % Compute color axes limits 
      
    CMin = min([real(V); imag(V)]);
    CMax = max([real(V); imag(V)]);
    CLim = [CMin CMax] + OFFSET*(CMax-CMin)*[-1 1];
    
    % Plot imaginary finite element solution  
      
    fig_1 = figure('Name','Crouzeix-Raviart finite elements');
    patch('faces', Elements, ...
          'vertices', [Coordinates(:,1) Coordinates(:,2) real(V)], ...
          'CData', real(V), ...
          'facecolor', 'interp', ...
          'edgecolor', 'none');
    set(gca,'XLim',XLim,'YLim',YLim,'CLim',CLim,'DataAspectRatio',[1 1 1]);
    fig_2 = figure('Name','Linear finite elements');
    patch('faces', Elements, ...
          'vertices', [Coordinates(:,1) Coordinates(:,2) imag(V)], ...
          'CData', imag(V), ...
          'facecolor', 'interp', ...
          'edgecolor', 'none');  
    set(gca,'XLim',XLim,'YLim',YLim,'CLim',CLim,'DataAspectRatio',[1 1 1]);
    
    if(nargout > 0)
      varargout{1} = fig_1;
      varargout{2} = fig_2;
    end
      
  end
  
return