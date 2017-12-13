function Aloc = assemDir_hpDG_1D(Coordinates,p,Ghandle,shap,grad_shap,s,alpha,varargin)
% ASSEMDIR_HPDG_1D assemble hpDG-FEM contributions.
%
%   A = ASSEMDir_HPDG_1D(COORDINATES,P,EHANDLE)
%
%   [I,J,A] = ASSEMDir_Bnd_HPDG(COORDINATES,P,EHANDLE)
%
%   [I,J,A] = ASSEMDir_Bnd_HPDG_1D(COORDINATES,P,EHANDLE,VARARGIN)
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
  
  Aloc=zeros(sum(p+1),1);
  % Preallocate memory
  
  nDofs = p(1)+p(N);
  I = zeros(nDofs,1);
  J = zeros(nDofs,1);
  A = zeros(nDofs,1);
  
  %left boundary
  
  %shap-functions
  g=Ghandle(Coordinates(1),varargin{:});
  h=Coordinates(2)-Coordinates(1);
  pm=p(1);
  Leg=shap(1,1:p(1)+1);
  GradLeg=2/h*grad_shap(1,1:p(1)+1);
  
  % Compute element stiffness matrix
    
  Aloc(1:p(1)+1) = s*GradLeg'*g+alpha*pm^2/h*Leg'*g;
  
  %right boundary
  
  %shap-functions
  g=Ghandle(Coordinates(N+1),varargin{:});
  h=Coordinates(N+1)-Coordinates(N);
  pm=p(N);
  Leg=shap(2,1:p(N)+1);
  GradLeg=2/h*grad_shap(2,1:p(N)+1);
  
  % Compute element stiffness matrix
    
  Aloc(end-p(N):end) = -s*GradLeg'*g+alpha*pm^2/h*Leg'*g;
  
return