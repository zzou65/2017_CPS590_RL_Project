function uex = uex(x,varargin)

%      u(r,theta) = r^2/3*sin(2/3*theta)


  % Compute polar coordinates

  r = sqrt(x(:,1).^2 + x(:,2).^2);
  theta = atan2(x(:,2),x(:,1));
  Loc = theta < 0;
  theta(Loc) = theta(Loc) + 2*pi;
  
  % Compute the exact solution for H1 semi norm
  
  uex = (r.^(2/3)).*sin(2/3*theta);
  
return