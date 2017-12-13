function Aloc = STIMA_InnPen_hpDG_1D(x,pl,pr,shap,alpha)
% STIMA_INNPEN_DGLFE Element stiffness matrix for interior terms.
%
%   ALOC = STIMA_INNPEN_hpDG_1D(x,pl,pr,varargin) computes the
%   entries of the element stiffness matrix for the interior terms.
%   
%   pl and pr are polynomial degrees on left and right element
%
%   Example:
%
%   Aloc = STIMA_InnPenn_hpDG_1D(x,pl,pr,1,varargin)
%
%   See also grad_shap_DGLFE.

%   Copyright 2007-2007 Patrick Meury, Holger Heumann
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
  rLeg=shap(1,1:pr+1);
  
  Aloc(1:(pl+1),1:(pl+1))=lLeg'*lLeg;
  Aloc(1:(pl+1),(pl+2):(pl+pr+2))=-lLeg'*rLeg;
  Aloc((pl+2):(pl+pr+2),1:(pl+1))=-rLeg'*lLeg;
  Aloc((pl+2):(pl+pr+2),(pl+2):(pl+pr+2))=rLeg'*rLeg;
  
 Aloc=alpha*max([pl,pr])^2/max([hl,hr])*Aloc;
 %Aloc=alpha*Aloc;
  
return