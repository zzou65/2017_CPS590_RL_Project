function varargout = plot_USR(Mesh)
% PLOT_USR Universal similarity region plot.
%
%   PLOT_USR(MESH) generates a plot of the universal similarity region of
%   all elements in the struct MESH.
%
%   H = PLOT_USR(MESH) also returns the handle to the generated figure.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh, where M
%                is equal to the number of vertices contained in the mesh.
%    ELEMENTS    M-by-3 matrix specifying the elements of the mesh, where M
%                is equal to the number of elements contained in the mesh. 
%   
%   Example:
%
%   plot_USR(Mesh);
%

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  nElements = size(Mesh.Elements,1);
  
  % Initialize constants
  
  MARKER_COLOR = 'k';
  MARKER_STYLE = 'x';
  USR_COLOR = 'k';
  USR_STYLE = '-';
  MID_STYLE = '-';
  STEPSIZE = 1e-4;
  XLIM = [-0.05 1.05];
  YLIM = [-0.05 1.05*sqrt(3)/2];
  
  % Compute USR coordinates
  
  X_USR = zeros(1,nElements);
  Y_USR = zeros(1,nElements);
  obtuse = 0;
  acute = 0;
  for i = 1:nElements
  
    % Extract current element from mesh  
      
    Coordinates = Mesh.Coordinates(Mesh.Elements(i,:),:);
    
    % Compute longest edge
    
    [h_max,id] = max([norm(Coordinates(3,:)-Coordinates(2,:)), ...
                      norm(Coordinates(1,:)-Coordinates(3,:)), ... 
                      norm(Coordinates(2,:)-Coordinates(1,:))]);
    
    % Translate back to origin and rescale the element
    
    j1 = rem(id+3,3)+1;
    X0 = Coordinates(j1,:);
    Coordinates = (Coordinates-repmat(X0,3,1))./h_max;
    
    % Rotate triangle and flip opposite vertex
    
    j2 = rem(id+4,3)+1; 
    theta = atan2(Coordinates(j2,2),Coordinates(j2,1));
    if(theta < 0)
      theta = theta+2*pi;
    end
    A = [cos(theta) -sin(theta); sin(theta) cos(theta)];
    Coordinates = Coordinates*A;
    if(Coordinates(id,2) < 0)
      X_USR(i) = Coordinates(id,1);
      Y_USR(i) = -Coordinates(id,2);
    else
      X_USR(i) = Coordinates(id,1);
      Y_USR(i) = Coordinates(id,2);
    end
    
    % Check wheter triangle is acute or obtuse
    
    if((X_USR(i)-0.5)^2+Y_USR(i)^2 > 1/4)
      acute = acute+1;  
    else
      obtuse = obtuse+1;  
    end
  end
  
  % Compute USR boundaries
  
  X_BD = 0:STEPSIZE:0.5;
  nSteps = size(X_BD,2);
  Y_MID = sqrt(1/4-(X_BD-1/2).^2);
  Y_UP = [sqrt(1-(X_BD-1).^2)];
  loc_1 = 2:nSteps;
  loc_2 = (nSteps-1):-1:1;
  X_BD = [X_BD X_BD(loc_1)+1/2];
  Y_MID = [Y_MID Y_MID(loc_2)];
  Y_UP = [Y_UP Y_UP(loc_2)];
  Y_LOW = zeros(1,2*nSteps-1);
    
  % Generate figure
  
  fig = figure('Name','Universal similarity region');
  hold on;
  plot(X_USR,Y_USR,[MARKER_COLOR MARKER_STYLE],'MarkerSize',3);
  plot(X_BD,Y_LOW,[USR_COLOR USR_STYLE], ...
       X_BD,Y_MID,[USR_COLOR MID_STYLE], ...
       X_BD,Y_UP,[USR_COLOR USR_STYLE]);
  hold off;
  title('{\bf Universal similarity region}');
  xlabel(['{\bf # Acute triangles  :  ',int2str(acute), ...
          ',      # Obtuse triangles  :  ',int2str(obtuse),'}'])
  box('on');
  set(gca,'XLim',XLIM,'YLim',YLIM);
  axis('equal');
  drawnow;
  
  % Assign output arguments
  
  if(nargout > 0)
    varargout{1} = fig;  
  end
  
return