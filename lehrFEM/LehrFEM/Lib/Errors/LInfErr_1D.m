function err = LInfErr_1D(Coordinates,u,FHandle,varargin)
% LINFERR_1D Discretization error in Linf norm for 1D linear finite
%            elements.
%
%   ERR = LINFERR_1D(COORDINATES,U,FHANDLE) computes the value of the
%   discretization error between the exact solution given by the function
%   handle FHANDLE and the finite element solution U.
%
%   ERR = LINFERR_1D(COORDINATES,U,FHANDLE,FPARAM) also handles the
%   variable length argument list FPARAM to the exact solution FHANDLE.
%
%   Example:
%
%   err = LInfErr_1D(Coordinates,u,FHandle);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
   
  err = max(abs(u-FHandle(Coordinates,varargin{:})));

return