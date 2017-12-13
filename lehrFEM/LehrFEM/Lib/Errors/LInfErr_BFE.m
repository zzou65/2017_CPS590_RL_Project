function err = LInfErr_BFE(Mesh,u,FHandle,varargin)
% LINFERR_BFE Discretization error in Linf norm for bilinear finite
%             elements.
%
%   ERR = LINFERR_BFE(MESH,U,FHANDLE) computes the discretization error
%   between the exact solution given by the function handle FHANDLE
%   and the finite element solution U on the struct MESH.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-4 matrix specifying the elements of the mesh. 
%
%   ERR = LINFERR_BFE(MESH,U,FHANDLE,FPARAM) also handles the variable
%   length argument list FPARAM to the exact solution FHANDLE.
%
%   Example:
%
%   err = LInfErr_BFE(Mesh,u,FHandle);

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
   
  err = max(abs(u-FHandle(Mesh.Coordinates,varargin{:})));

return