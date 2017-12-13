function U = L2_interp_QFE(Mesh,Quadrule,F_Handle,varargin)
% L2_INTERP_QFE L2 projection with QFE
%
% L2_INTERP_QFE(MESH,QUADRULE,F_HANDLE) computes the solution of the L2 
% projection using quadratic finite element method.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%    EDGES        P-by-2 matrix specifying all edges of the mesh.
%    VERT2EDGE    M-by-M sparse matrix which specifies wheter the two
%                 vertices i and j are connected by an edge with number
%                 VERT2EDGE(i,j).
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used to
%   do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%
%   Example:
%
%   U = L2_interp_QFE(Mesh,P7O6(),F_Handle);
%
%   See also L2_interp_LFE.

    % Assemble stiffness matrix and mass matrix
    
    M = assemMat_QFE(Mesh,@MASS_QFE,P7O6());
    L = assemLoad_QFE(Mesh,Quadrule,F_Handle);
    
    % Solve the solution
    
    U = M\L;
    
return    
    
