function uint = comp_interp(Mesh,CDofs,EDofs,QuadRule,Shap,UHandle,varargin)
% computes the L^2 projection of the solution to the finite element space
% for hp finite elements

% Copyright 2009 Christoph Wiesmeyr
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland


Elem2Dof = build_DofMaps(Mesh,EDofs,CDofs);

% compute the right hand side
L = assemLoad_hp(Mesh,Elem2Dof,QuadRule,Shap,UHandle);

% compute the mass matrix
M = assemMat_hp(Mesh,Elem2Dof,@MASS_hp,QuadRule,Shap);

% solve linear system
uint = M\L;