function Aloc = STIMA_Inn_Lapl_hpDG_1D(x,pl,pr,shap,grad_shap,s)
% STIMA_INN_Lapl_DGLFE Element stiffness matrix for interior terms.
%
%   ALOC = STIMA_INN_Lapl_hpDG_1D(x,pl,pr,varargin) computes the
%   entries of the element stiffness matrix for the interior terms.
%   
%   pl and pr are polynomial degrees on left and right element
%
%   The integer S can specifies wheter the diffusive fluxes are discretized
%   in a symmetric or anti-symmetric way:
%    +1 Antisymmetric discretization of diffusive fluxes
%    -1 Symmetric discretization of diffusive fluxes
%
%   Example:
%
%   Aloc = STIMA_Inn_Lapl_hpDG_1D(x,pl,pr,1,varargin)
%
%   See also grad_shap_DGLFE.

%   Copyright 2007-2007 Patrick Meury Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Preallocate memory
  
  Aloc = zeros(pl+pr+2,pl+pr+2);

  % left and right elementsize
  hl=x(2)-x(1);
  hr=x(3)-x(2);
  
  % values of basis functions and their gradients on element boundary
  lLeg=shap(2,1:pl+1);
  lGradLeg=2/hl*grad_shap(2,1:pl+1);
  rLeg=shap(1,1:pr+1);
  rGradLeg=2/hr*grad_shap(1,1:pr+1);
  
  Aloc(1:(pl+1),1:(pl+1))=lLeg'*lGradLeg+s*lGradLeg'*lLeg;
  Aloc((pl+2):(pl+pr+2),1:(pl+1))=-rLeg'*lGradLeg+s*rGradLeg'*lLeg;
  Aloc(1:(pl+1),(pl+2):(pl+pr+2))=+lLeg'*rGradLeg-s*lGradLeg'*rLeg;
  Aloc((pl+2):(pl+pr+2),(pl+2):(pl+pr+2))=-rLeg'*rGradLeg-s*rGradLeg'*rLeg;
  
  Aloc=1/2*(Aloc);
  
return