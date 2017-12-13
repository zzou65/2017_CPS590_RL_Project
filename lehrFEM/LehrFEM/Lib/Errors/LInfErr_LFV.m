function err = LInfErr_LFV(Mesh,u,FHandle,varargin)
% LINFERR_LFV Discretization error in Linf norm for linear finite volumes.
%
%   ERR = LINFERR_LFV(MESH,U,FHANDLE) computes the discretization error
%   between the exact solution given by the function handle FHANDLE
%   and the finite volume solution U on the struct MESH.
%
%   The struct MESH should at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 matrix specifying the elements of the mesh. 
%
%   ERR = LINFERR_LFV(MESH,U,FHANDLE,FPARAM) also handles the variable
%   length argument list FPARAM to the exact solution FHANDLE.
%
%   Example:
%
%   err = LInfErr_LFV(Mesh,u,fhandle);

%   Copyright 2007-2007 Eivind Fonn
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
   
    err = 0;
    for i=1:size(Mesh.Coordinates,1)
        err = max(err, abs(u(i)-FHandle(Mesh.Coordinates(i,:))));
    end

return
