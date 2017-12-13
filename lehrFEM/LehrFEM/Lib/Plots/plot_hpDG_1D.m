function varargout = plot_hpDG_1D(Coordinates,p,u,Shap)
% PLOT_HPDG_1D Plot finite element solution.
%
%   PLOT_HPDG_1D(CORDINATES,P,U,SHAP) generates a plot of the hpDG solution
%   U on the mesh specified by COORDINATES using the polynomial degrees
%   specified by P and the shape functions provided by the function handle
%   SHAP on each element.
%
%   H = PLOT_HPDG_1D(COORDINATES,P,U,SHAP) also returns the handle to the
%   figure.
%
%   Example:
%
%   PLOT_hpDG(coordinates, p, u, Leg);

%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  OFFSET = 0.05;
  
  % Compute axes limits
  
  nCoordinates = size(Coordinates,2);
  nElements = nCoordinates-1;
  
  XMin = min(Coordinates);
  XMax = max(Coordinates);
  XLim = [XMin XMax] + OFFSET*(XMax-XMin)*[-1 1]; 
    
  % Evaluate shape functions on reference element
  
  pmax = max(p);
  nPts = 2*(pmax+1);
  xhat = transpose(-1 + 2*(0:nPts)/nPts);
  yhat = Shap(xhat,pmax);
    
  % Evaluate hpDG solution on each element
  
  x = zeros(2,nElements*nPts);
  y = zeros(2,nElements*nPts);
 
  offset = 0;
  for i=1:nElements
    
    % Compute element mapping  
      
    xx = (Coordinates(i+1)+Coordinates(i))/2 + ...
         (Coordinates(i+1)-Coordinates(i))/2*xhat;
    
    % Evaluate solution
    
    nDofs = p(i)+1;
    yy = zeros(size(xhat));
    for j = 1:nDofs
      yy = yy + u(offset+j)*yhat(:,j);
    end
    offset = offset + nDofs;
    
    x_start(1,:) = xx(1:(nPts));
    x_start(2,:) = xx(2:(nPts+1));
    x_end(1,:) = yy(1:nPts);
    x_end(2,:) = yy(2:(nPts+1));
    
    loc = (i-1)*nPts + (1:nPts);
    x(:,loc) = x_start;
    y(:,loc) = x_end;
  
  end
  
  % Compute axes limits
  
  YMin = min(min(y));
  YMax = max(max(y));
  if(YMin < YMax)
    YLim = [YMin YMax] + OFFSET*(YMax-YMin)*[-1 1];
  else
    YLim = [1-OFFSET 1+OFFSET]*YMin;   
  end

  % Generate figure
  
  fig = figure('Name','Linear finite elements');
  plot(x,y,'r-');
  set(gca,'XLim',XLim,'YLim',YLim);
    
  % Assign output arguments
  
  if(nargout > 0)
    varargout{1} = fig;  
  end
  
return