function varargout = plot_1D(U,Coordinates)
% PLOT_1D Plot finite element solution.
%
%   PLOT_1D(U,COORDINATES) generates a plot for the finite element solution
%   U on the mesh COORDINATES.
%
%   H = PLOT_1D(U,MESH) also returns the handle to the figure.
%
%   Example:
%
%   plot_1D(U,Coordinates);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  OFFSET = 0.05;
  
  % Compute axes limits
  
  XMin = min(Coordinates);
  XMax = max(Coordinates);
  XLim = [XMin XMax] + OFFSET*(XMax-XMin)*[-1 1]; 
    
  % Generate figure
    
  if(isreal(U))
  
    % Compute color axes limits 
      
    YMin = min(U);
    YMax = max(U);
    if(YMin < YMax)
      YLim = [YMin YMax] + OFFSET*(YMax-YMin)*[-1 1];
    else
      YLim = [1-OFFSET 1+OFFSET]*YMin;   
    end
         
    % Plot real finite element solution  
      
    fig = figure('Name','Linear finite elements');
    plot(Coordinates,U,'r-');
    set(gca,'XLim',XLim,'YLim',YLim,'DataAspectRatio',[1 1 1]);
    
    if(nargout > 0)
      varargout{1} = fig;
    end
 
  else  
      
    % Compute color axes limits 
      
    YMin = min([real(U); imag(U)]);
    YMax = max([real(U); imag(U)]);
    if(YMin < YMax)
      YLim = [YMin YMax] + OFFSET*(YMax-YMin)*[-1 1];
    else
      YLim = [1-OFFSET 1+OFFSET]*YMin;   
    end
    
    % Plot imaginary finite element solution  
      
    fig_1 = figure('Name','Linear finite elements');
    plot(Coordinates,real(U),'r-');
    set(gca,'XLim',XLim,'YLim',YLim,'DataAspectRatio',[1 1 1]);
    fig_2 = figure('Name','Linear finite elements');
    plot(Coordinates,imag(U),'r-');
    set(gca,'XLim',XLim,'YLim',YLim,'DataAspectRatio',[1 1 1]);
    if(nargout > 0)
      varargout{1} = fig_1;
      varargout{2} = fig_2;
    end
      
  end
  
return