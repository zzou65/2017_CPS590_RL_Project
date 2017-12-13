function varargout = assemMat_Vol_hpDG_1D(Coordinates,p,EHandle,varargin)
% ASSEMMAT_VOL_HPDG_1D Assemble hpDG-FEM contributions.
%
%   A = ASSEMMAT_VOL_HPDG_1D(COORDINATES,P,EHANDLE)
%
%   [I,J,A] = ASSEMMAT_VOL_HPDG(COORDINATES,P,EHANDLE)
%
%   [I,J,A] = ASSEMMAT_VOL_HPDG_1D(COORDINATES,P,EHANDLE,VARARGIN)

%   Copyright 2007-2007 Patrick Meury & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nCoordinates = size(Coordinates,2);
  nElements = nCoordinates-1;
  
  % Preallocate memory
  
  nDofs = sum((p+1).^2);
  I = zeros(nDofs,1);
  J = zeros(nDofs,1);
  A = zeros(nDofs,1);
  
  % Assemble matrix
  
  offset_idx = 0; 
  offset_dof = 0;
  vidx = [1 2];
  for i = 1:nElements
      
    % Compute element stiffness matrix
    Aloc = EHandle(Coordinates(vidx),p(i),varargin{:});
      
    % Add contributions to global vector
      
    nDofs = p(i)+1;
    idx = offset_dof + (1:nDofs);
    loc = offset_idx + (1:(nDofs^2));
    
    I(loc) = set_Rows(idx,nDofs);
    J(loc) = set_Cols(idx,nDofs);
    A(loc) = Aloc(:);
    
    % Update counters
    
    vidx = vidx+1;
    offset_dof = offset_dof + nDofs;
    offset_idx = offset_idx + nDofs^2;
    
  end
  
  % Assign output arguments
  
  if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
  else
    varargout{1} = sparse(I,J,A);
  end
  
return