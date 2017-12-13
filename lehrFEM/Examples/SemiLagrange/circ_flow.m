function flow = circ_flow(x,varargin)

%   Copyright 2008-2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants

  R0 = 0.9;       % Radius
  C = [0 0];  % Center
  
  % Preallocate memory
  
  flow = zeros(size(x,1),2);
    
  r = sqrt((x(:,1)-C(1)).^2+(x(:,2)-C(2)).^2);
  Loc = find(r < R0);
  flow(Loc,:) = [-x(Loc,2) +x(Loc,1)];
  
return