function err = TVDNorm_hpDG_1D(Coordinates,p,u,QuadRule,Shap,varargin)
% TVD-norm for 1D hpDG
% 
%   Coordinates: n+1 Meshpoints
%             u: solution for last step
%             p: polynomial approxiamtion each of the n cells
%          shap: basisfunctions evaluated at quadrature points
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
  
  % Compute L2 error
  
  err = 0;
  offset = 0;
  u_aver=zeros(nElements,1);
  for i = 1:nElements
     
    % Compute finite element solution
      
    u_h = zeros(size(QuadRule.x));
    nDofs = p(i)+1;
    for j = 1:nDofs
      u_h = u_h + u(offset+j)*Shap(:,j);  
    end
    
    % Compute element map
    
    h = (Coordinates(i+1)-Coordinates(i))/2;
     
    % Compute L2 error on current element
    
    u_aver(i)=sum(QuadRule.w.*u_h)*h;
    
    % Update counter
    
    offset = offset + nDofs;
    
  end
  
  err = sum(abs(u_aver(2:end)-u_aver(1:(end-1))));

return