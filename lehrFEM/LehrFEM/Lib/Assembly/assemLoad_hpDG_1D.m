function L = assemLoad_hpDG_1D(Coordinates,p,QuadRule,Shap,FHandle,varargin)
% ASSEMLOAD_HPDG_1D Assemble hpDG-FEM contributions.
%
%   L = ASSEMLOAD_HPDG_1D(COORDINATES,P,QUADRULE,SHAP,FHANDLE)
%
%   L = ASSEMLOAD_HPDG_1D(COORDINATES,P,QUADRULE,SHAP,FHANDLE,VARARGIN)
%

%   Copyright 2007-2007 Patrick Meury & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nCoordinates = size(Coordinates,2);
  nElements = nCoordinates-1;

  % Preallocate memory
  
  nDofs = sum(p)+nElements;
  L = zeros(nDofs,1);
  
  % Assemble load vector
  
  offset = 0;
  for i = 1:nElements
  
    % Evaluate element mapping
    
    h = (Coordinates(i+1)-Coordinates(i))/2;
    x = (Coordinates(i+1)+Coordinates(i))/2 + h*QuadRule.x;
    
    % Evaluate load data
    
    F = FHandle(x,varargin{:});
      
    % Evaluate element load vector
    
    nDofs = p(i)+1;
    for i = 1:nDofs
      L(offset+i) = L(offset+i) + sum(QuadRule.w.*F.*Shap(:,i))*h;
    end
    
    % Update counter
    
    offset = offset+nDofs;    
    
  end
  
return