function U = LFE_DOF_interp(Mesh,F_Handle,varargin)
% LFE_DOF_interp linear interpolation using LFE DOFs
%
% LFE_DOF_interp(MESH,F_HANDLE) returns the LFe DOFS of function in F_HANDLE
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%
%   
%   Example:
%
%   U = LFE_DOF_interp(Mesh,F_Handle);
%
%   See also QFE_DOF_interp.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    U = F_Handle(Mesh.Coordinates);
    
return    
    
