function U = H1_interp_LFE(Mesh,Quadrule,F_Handle,GD_Handle,varargin)
% H1_INTERP_LFE H1 projection with LFE
%
% H1_INTERP_LFE(MESH,QUADRULE,F_HANDLE,GD_HANDLE) computes the solution of 
% the H1 projection using linear finite element method.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used to
%   do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%
%   Example:
%
%   U = H1_interp_LFE(Mesh,P7O6(),F_Handle,GD_Handle);
%
%   See also H1_interp_QFE.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Assemble stiffness matrix and mass matrix
    
    [IA JA KA] = assemMat_LFE(Mesh,@STIMA_Lapl_LFE,P7O6());
    [IM JM KM] = assemMat_LFE(Mesh,@MASS_LFE,P7O6());
    L = assemInterp_LFE(Mesh,Quadrule,F_Handle,GD_Handle);
    
    % Solve the solution
    
    U = sparse([IA ; IM],[JA ; JM],[KA ; KM])\L;
    
return    
    
