function varargout = assemMat_Inn_hpDG_1D(Coordinates,p,EHandle,varargin)
% ASSEMMAT_Inn_HPDG_1D assemble hpDG-FEM contributions.
%
%   A = ASSEMMAT_Inn_HPDG_1D(COORDINATES,P,EHANDLE)
%
%   [I,J,A] = ASSEMMAT_Inn_HPDG(COORDINATES,P,EHANDLE)
%
%   [I,J,A] = ASSEMMAT_Inn_HPDG_1D(COORDINATES,P,EHANDLE,VARARGIN)
%
%    p: n-1 array with polynomial degree of each element
%    Coordinates: n-1 array specifying the Gridpoints
%    EHandle: function handle for Contribution 

%   Copyright 2007-2007 Patrick Meury Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nCoordinates = size(Coordinates,2);
  nElements = nCoordinates-1;
  
  % Preallocate memory
  
  nDofs = sum((p+1));
  I = zeros(nDofs,1);
  J = zeros(nDofs,1);
  A = zeros(nDofs,1);
  
  % Assemble matrix
  
  offset_idx = 0; 
  offset_dof = 0;
  for i = 2:(nCoordinates-1)
      
    % Compute element stiffness matrix
    
    Aloc = EHandle(Coordinates(i-1:i+1),p(i-1),p(i),varargin{:});
      
    % Add contributions to global vector
    
    ln_idx = p(i-1)+1;
    rn_idx = p(i)+1;
    n_idx=ln_idx+rn_idx;
    idx = offset_dof + [1:n_idx];
    
    loc = offset_idx + [1:n_idx^2];
    I(loc) = set_Rows(idx,n_idx);
    J(loc) = set_Cols(idx,n_idx);
    A(loc) = Aloc(:);
    
    % Update counters
    
    offset_dof = offset_dof + ln_idx;
    offset_idx = offset_idx + n_idx^2;
    
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