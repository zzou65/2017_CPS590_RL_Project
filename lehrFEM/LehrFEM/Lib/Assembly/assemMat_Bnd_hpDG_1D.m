function varargout = assemMat_Bnd_hpDG_1D(Coordinates,p,shap,grad_shap,s,alpha)
% ASSEMMAT_Bnd_HPDG_1D assemble hpDG-FEM contributions.
%
%   A = ASSEMMAT_Bnd_HPDG_1D(COORDINATES,P,EHANDLE)
%
%   [I,J,A] = ASSEMMAT_Bnd_HPDG(COORDINATES,P,EHANDLE)
%
%   [I,J,A] = ASSEMMAT_Bnd_HPDG_1D(COORDINATES,P,EHANDLE,VARARGIN)
%
%    p: n-1 array with polynomial degree of each element
%    Coordinates: n-1 array specifying the Gridpoints
%    EHandle: function handle for Contribution 

%   Copyright 2007-2007 Patrick Meury Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  N=size(Coordinates,2)-1;
  offset_dof=0;
  offset_idx=0;
  % Preallocate memory
  
  nDofs = (p(1)+1)^2+(p(N)+1)^2;
  I = zeros(nDofs,1);
  J = zeros(nDofs,1);
  A = zeros(nDofs,1);
  
  %left boundary
  
  %shap-functions
  h=Coordinates(2)-Coordinates(1);
  pm=p(1);
  Leg=shap(1,1:p(1)+1);
  GradLeg=2/h*grad_shap(1,1:p(1)+1);
  
  % Compute element stiffness matrix
    
  Aloc = -Leg'*GradLeg-s*GradLeg'*Leg-alpha*pm^2/h*Leg'*Leg;
  
  % Add contributions to global vector
    
  n_idx=p(1)+1;
  idx = [1:n_idx];
    
  loc = [1:n_idx^2];
  I(loc) = set_Rows(idx,n_idx);
  J(loc) = set_Cols(idx,n_idx);
  A(loc) = Aloc(:);
    
  %right boundary
  
  %shap-function
  h=Coordinates(N+1)-Coordinates(N);
  pm=p(N);
  Leg=shap(2,1:p(N)+1);
  GradLeg=2/h*grad_shap(2,1:p(N)+1);
  
  % Compute element stiffness matrix
    
  Aloc = (+Leg'*GradLeg+s*GradLeg'*Leg)-alpha*pm^2/h*Leg'*Leg;
  
  % Add contributions to global vector
    
  offset=(p(1)+1)^2;
  idx =sum(p+1)-p(N):sum(p+1);
  n_idx=(p(N)+1);
  
  loc =offset+[1:n_idx^2];
  I(loc) = set_Rows(idx,n_idx);
  J(loc) = set_Cols(idx,n_idx);
  A(loc) = Aloc(:);
    
  % Assign output arguments
  
  if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
  else
    varargout{1} = sparse(I,J,A,sum(p+1),sum(p+1));
  end
  
return