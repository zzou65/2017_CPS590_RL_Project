function err = LInfErr_QFE(Mesh,u,FHandle,varargin)
% LINFERR_QFE Discretization error in Linf norm for quadratic finite
%             elements.
%
%   ERR = LINFERR_QFE(MESH,U,FHANDLE) computes the discretization error
%   between the exact solution given by the function handle FHANDLE and the
%   finite element solution U on the struct MESH.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 matrix specifying the elements of the mesh. 
%    EDGES       P-by-2 matrix specifying the edges of the mesh.
%
%   ERR = LINFERR_QFE(MESH,U,FHANDLE,FPARAM) also handles the variable
%   length argument list FPARAM to the exact solution FHANDLE.
%
%   Example:
%
%   err = LInfErr_QFE(Mesh,u,fhandle);

%   Copyright 2005-2005 Patrick Meury & Kah Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
   
  x_mid = 1/2*(Mesh.Coordinates(Mesh.Edges(:,1),:) + ...
               Mesh.Coordinates(Mesh.Edges(:,2),:));
  err = max(abs(u-FHandle([Mesh.Coordinates; x_mid],varargin{:})));
               
return
