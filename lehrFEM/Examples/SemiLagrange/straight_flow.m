function flow = straight_flow(x,C,R0,varargin)

%   Copyright 2008-2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants

  %R0 = 0.4;       % Radius
  %C = [0.5 0.5];  % Center
  
  % Preallocate memory
  
  flow = zeros(size(x,1),2);
    
  r = max(abs(x-ones(size(x,1),1)*C),[],2);%sqrt((x(:,1)-C(1)).^2+(x(:,2)-C(2)).^2);
  Loc = find(r < R0);
  flow(Loc,:) = sqrt(2)*[ones(size(Loc,1),1) 1*ones(size(Loc,1),1)];
  
return