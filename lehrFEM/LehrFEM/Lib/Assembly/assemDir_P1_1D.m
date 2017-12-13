function [U,FreeDofs] = assemDir_P1_1D(Coordinates,DNodes,FHandle,varargin)
% ASSEMDIR_1D Dirichlet boundary conditions.
%
%   [U,FREEDOFS] = ASSEMDIR_P1_1D(COORDINATES,DNODES,FHANDLE) incoporates
%   the Dirichlet boundary conditions with the data given by FHANDLE at the
%   vertices denoted by DNODES into the finite element solution U.
%
%   FREEDOFS is a matrix specifying the vertices with no prescribed
%   Dirichlet boundary data.  
%
%   [U,FREEDOFS] = ASSEMDIR_P1_1D(COORDINATES,DNODES,FHANDLE,FPARAM) also
%   handles the variable length argument list FPARAM to the boundary data
%   function FHANDLE.
%
%   Example:
%
%   [U,FreeDofs] = assemDir_P1_1D(Coordiantes,[1 size(Coordinates,1)],GD_Handle);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants

  nCoordinates = size(Coordinates,1);
  
  % Compute set of free dofs 
 
  FreeDofs = setdiff(1:nCoordinates,DNodes);
  
  % Incorporate Dirichlet boundary condition
  
  U = zeros(nCoordinates,1);
  U(DNodes(1)) = FHandle(Coordinates(DNodes(1)),DNodes(1),varargin{:});
  if(size(DNodes,1) > 1)
    U(DNodes(2)) = FHandle(Coordinates(DNodes(2)),DNodes(2),varargin{:});    
  end
   
return