function U = L2_interp_LFE(Mesh,Quadrule,F_Handle,varargin)
% L2_INTERP_LFE L2 projection with LFE
%
% L2_INTERP_LFE(MESH,QUADRULE,F_HANDLE) computes the solution of the L2 
% projection using linear finite element method.
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
%   U = L2_interp_LFE(Mesh,P7O6(),F_Handle);
%
%   See also L2_interp_QFE.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Assemble stiffness matrix and mass matrix
    
    M = assemMat_LFE(Mesh,@MASS_LFE,P7O6());
    L = assemLoad_LFE(Mesh,Quadrule,F_Handle);
    
    % Solve the solution
    
    U = M\L;
    
return    
    
