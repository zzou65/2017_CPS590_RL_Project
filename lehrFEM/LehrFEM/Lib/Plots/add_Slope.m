function [] = add_Slope(axes_handle,loc,rate,varargin)
% ADD_SLOPE Add slope triangle to current plot.
%
%   ADD_SLOPE(AXES_HANDLE,LOC,RATE) adds a slope triangle with convergence
%   rate P at the location LOC to the plot specified by AXES_HANDLE using
%   color specified in varargin, default is black.
%
%   The string LOC specifies the placement of the slope triangle:
%    'North'      Placed at top.
%    'NorthEast'  Placed into upper right-hand corner.
%    'East'       Placed on right-hand side.
%    'SouthEast'  Placed into lower right-hand corner.
%    'South'      Placed at bottom.
%    'SouthWest'  Placed into lower left-hand corner.
%    'West'       Placed on left-hand side.
%    'NorthWest'  Placed into upper left-hand corner.
%
%   Example:
%
%   add_Slope(gca,'NorthEast',1/2);
  
%   Copyright 2005-2006 Patrick Meury & Ralf Hiptmair & Holger Heumann
%   SAM -  Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  F_SIZE = 10;  % Font size
  
  % Get axes limits
  
  XLim = get(axes_handle,'XLim');
  YLim = get(axes_handle,'YLim');

  % Compute triangle location
  
  if(strcmp(lower(loc),'north') || strcmp(lower(loc),'n'))      
    
    % Compute relative coordinates  
      
    rl = 0.4;
    rr = 0.6;
    ry = 0.9;
    
    % Compute vertices of triangle
    
    xr = (XLim(2)/XLim(1))^rr*XLim(1);
    xl = (XLim(2)/XLim(1))^rl*XLim(1);
    xmid = sqrt(xl*xr);
    if(rate > 0)
      yr = (YLim(2)/YLim(1))^ry*YLim(1);
      yl = yr*(xl/xr)^rate;
      valign = 'top';
    else
      yl = (YLim(2)/YLim(1))^ry*YLim(1);
      yr = yl*(xr/xl)^rate;
      valign = 'bottom';
    end
    ymid = yl;
    
  elseif(strcmp(lower(loc),'northeast') || strcmp(lower(loc),'ne'))
    
    % Compute relative coordinates
      
    rl = 0.7;
    rr = 0.9;
    ry = 0.9;
    
    % Compute vertices of triangle
    
    xr = (XLim(2)/XLim(1))^rr*XLim(1);
    xl = (XLim(2)/XLim(1))^rl*XLim(1);
    xmid = sqrt(xl*xr);
    if(rate > 0)
      yr = (YLim(2)/YLim(1))^ry*YLim(1);
      yl = yr*(xl/xr)^rate;
      valign = 'top';
    else
      yl = (YLim(2)/YLim(1))^ry*YLim(1);
      yr = yl*(xr/xl)^rate;
      valign = 'bottom';
    end
    ymid = yl;
    
  elseif(strcmp(lower(loc),'east') || strcmp(lower(loc),'e'))
      
    % Compute relative coordinates
      
    rl = 0.7;
    rr = 0.9;
    ry = 0.6;
    
    % Compute vertices of triangle
    
    xr = (XLim(2)/XLim(1))^rr*XLim(1);
    xl = (XLim(2)/XLim(1))^rl*XLim(1);
    xmid = sqrt(xl*xr);
    if(rate > 0)
      yr = (YLim(2)/YLim(1))^ry*YLim(1);
      yl = yr*(xl/xr)^rate;
      valign = 'top';
    else
      yl = (YLim(2)/YLim(1))^ry*YLim(1);
      yr = yl*(xr/xl)^rate;
      valign = 'bottom';
    end
    ymid = yl;
    
  elseif(strcmp(lower(loc),'southeast') || strcmp(lower(loc),'se'))
    
    % Compute relative coordinates  
      
    rl = 0.7;
    rr = 0.9;
    ry = 0.1;
    
    % Compute vertices of triangle
    
    xr = (XLim(2)/XLim(1))^rr*XLim(1);
    xl = (XLim(2)/XLim(1))^rl*XLim(1);
    xmid = sqrt(xl*xr);
    if(rate > 0)
      yl = (YLim(2)/YLim(1))^ry*YLim(1);
      yr = yl*(xr/xl)^rate;
      valign = 'top';
    else
      yr = (YLim(2)/YLim(1))^ry*YLim(1);
      yl = yr*(xl/xr)^rate;
      valign = 'bottom';
    end
    ymid = yl;
    
  elseif(strcmp(lower(loc),'south') || strcmp(lower(loc),'s'))
    
    % Compute relative coordinates
      
    rl = 0.4;
    rr = 0.6;
    ry = 0.1;
    
    % Compute vertices of triangle
    
    xr = (XLim(2)/XLim(1))^rr*XLim(1);
    xl = (XLim(2)/XLim(1))^rl*XLim(1);
    xmid = sqrt(xl*xr);
    if(rate > 0)
      yl = (YLim(2)/YLim(1))^ry*YLim(1);
      yr = yl*(xr/xl)^rate;
      valign = 'top';
    else
      yr = (YLim(2)/YLim(1))^ry*YLim(1);
      yl = yr*(xl/xr)^rate;
      valign = 'bottom';
    end
    ymid = yl;
    
  elseif(strcmp(lower(loc),'southwest') || strcmp(lower(loc),'sw'))
    
    % Compute relative coordinates  
      
    rl = 0.1;
    rr = 0.3;
    ry = 0.1;
    
    % Compute vertices of triangle
    
    xr = (XLim(2)/XLim(1))^rr*XLim(1);
    xl = (XLim(2)/XLim(1))^rl*XLim(1);
    xmid = sqrt(xl*xr);
    if(rate > 0)
      yl = (YLim(2)/YLim(1))^ry*YLim(1);
      yr = yl*(xr/xl)^rate;
      valign = 'top';
    else
      yr = (YLim(2)/YLim(1))^ry*YLim(1);
      yl = yr*(xl/xr)^rate;
      valign = 'bottom';
    end
    ymid = yl;
    
  elseif(strcmp(lower(loc),'west') || strcmp(lower(loc),'w'))
    
    % Compute relative coordinates  
      
    rl = 0.1;
    rr = 0.3;
    ry = 0.6;
    
    % Compute vertices of triangle
    
    xr = (XLim(2)/XLim(1))^rr*XLim(1);
    xl = (XLim(2)/XLim(1))^rl*XLim(1);
    xmid = sqrt(xl*xr);
    if(rate > 0)
      yr = (YLim(2)/YLim(1))^ry*YLim(1);
      yl = yr*(xl/xr)^rate;
      valign = 'top';
    else
      yl = (YLim(2)/YLim(1))^ry*YLim(1);
      yr = yl*(xr/xl)^rate;
      valign = 'bottom';
    end
    ymid = yl;
    
  elseif(strcmp(lower(loc),'northwest') || strcmp(lower(loc),'nw'))
    
    % Compute relative coordinates  
      
    rl = 0.1;
    rr = 0.3;
    ry = 0.9;
    
    % Compute vertices of triangle
    
    xr = (XLim(2)/XLim(1))^rr*XLim(1);
    xl = (XLim(2)/XLim(1))^rl*XLim(1);
    xmid = sqrt(xl*xr);
    if(rate > 0)
      yr = (YLim(2)/YLim(1))^ry*YLim(1);
      yl = yr*(xl/xr)^rate;
      valign = 'top';
    else
      yl = (YLim(2)/YLim(1))^ry*YLim(1);
      yr = yl*(xr/xl)^rate;
      valign = 'bottom';
    end
    ymid = yl;
    
  else
    error('Invalid location.');  
  end
  
  % Add slope triangle to current plot
  
  if(~ishold)
    hold on;    
  end
  
  if (nargin==3)
   plot([xl xr xr xl],[yl yl yr yl],'k-');
   text(xmid,ymid,sprintf('p = %1.2f',rate), ...
       'FontSize',F_SIZE, ...
       'HorizontalAlignment','center', ...
       'VerticalAlignment',valign);
   hold off;
  end
  if (nargin==4)
   plot([xl xr xr xl],[yl yl yr yl],varargin{:});
   text(xmid,ymid,sprintf('p = %1.2f',rate), ...
       'FontSize',F_SIZE, ...
       'HorizontalAlignment','center', ...
       'VerticalAlignment',valign);
   hold off;
  end
return