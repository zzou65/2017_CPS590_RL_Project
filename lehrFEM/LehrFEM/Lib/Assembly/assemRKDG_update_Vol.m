function L = assemRKDG_update_Vol(Coordinates,u,p,shap,gradshap,flowhandle,QuadRule,varargin)
% RKDG_update_Inn
% Update for Volume term for RKDG timestepping
%
%   Coordinates: n+1 Meshpoints
%             u: solution for last step
%             p: polynomial approxiamtion each of the n cells
%          shap: basisfunctions evaluated at quadrature points
%      gradshap: gradient of basis functions evaluated at quadrature points
%    flowhandle: flux function
%      Quadrule: 1D quadrature rule 
%
%   Example:
%
%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

 % Initialize constants
  
  nCoordinates = size(Coordinates,2);
  nElements = nCoordinates-1;
  nDofs = sum(p)+nElements;
  
  % allocate memory
  
  L = zeros(nDofs,1);
  
  % Assemble load vector
  
  offset = 0;
  for i = 1:nElements
  
    % Evaluate element mapping
    
    h = (Coordinates(i+1)-Coordinates(i))/2;
    x = (Coordinates(i+1)+Coordinates(i))/2 + h*QuadRule.x;
    
    % evaluate u(x) and f(u)
    
    u_cell=shap(:,1:p(i)+1)*u(offset+(1:p(i)+1));
    F = flowhandle(u_cell,varargin{:});
      
    nDofs = p(i)+1;
    for j = 1:nDofs
      L(offset+j) = L(offset+j) + sum(QuadRule.w.*F.*gradshap(:,j));
    end
    
    % Update counter
    
    offset = offset+nDofs;    
    
  end
  
return