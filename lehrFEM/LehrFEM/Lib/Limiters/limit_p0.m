function u_new = limit_p0(Mesh,u_old,vargin)
% LIMIT_P0 Non-adaptive, piecewise constant limiter.
%
%   U = LIMIT_P0(MESH,U) projects the piecewie linear function U onto the
%   space of piecewise constant functions.
%
%   Example:
%
%   U = limit_p0(Mesh,U);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nElements = size(Mesh.Elements,1);  % Number of elements
  
  % Prealloctae memory
  
  u_new = zeros(3*nElements,1);
  
  % Compute limited value
  
  idx = [1 2 3];
  for i = 1:nElements
    u_new(idx) = sum(u_old(idx))/3;
    idx = idx+3;
  end
  
return