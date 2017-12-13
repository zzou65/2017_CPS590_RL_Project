function err = LInfErr_LFE(Mesh,u,FHandle,varargin)
% LINFERR_LFE Discretization error in Linf norm for linear finite elements.
%
%   ERR = LINFERR_LFE(MESH,U,FHANDLE) computes the discretization error
%   between the exact solution given by the function handle FHANDLE
%   and the finite element solution U on the struct MESH.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 matrix specifying the elements of the mesh. 
%
%   ERR = LINFERR_LFE(MESH,U,FHANDLE,FPARAM) also handles the variable
%   length argument list FPARAM to the exact solution FHANDLE.
%
%   Example:
%
%   err = LInfErr_LFE(Mesh,u,fhandle);

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
   
  err = max(abs(u-FHandle(Mesh.Coordinates,varargin{:})));

return
