function L = assemNeu_P1_1D(Coordinates,NNodes,L,FHandle,varargin)
% ASSEMDIR_1D Neumann boundary conditions.
%
%   L = ASSEMNEU_P1_1D(COORDINATES,NNODES,L,FHANDLE) adds the Neumann
%   boundary conditions with the data given by FHANDLE at the vertices
%   denoted by NNODES onto the right-hand side load vector L and the
%   finite element solution U.
%
%   L = ASSEMNEU_P1_1D(COORDINATES,NNODES,L,FHANDLE,FPARAM) also handles
%   the variable length argument list FPARAM to the boundary data function
%   FHANDLE.
%
%   Example:
%
%   Coordinates = transpose(0:.01:1);
%   F_Handle = @(x,varargin)8*ones(size(x));
%   GN_Handle = @(x,varargin)zeros(size(x));
%   A = assemMat_P1_1D(Coordinates,@STIMA_Lapl_P1_1D);
%   L = assemLoad_P1_1D(Coordinates,@LOAD_1D,gauleg(0,1,5),F_Handle);
%   NNodes = 1;
%   L = assemNeu_P1_1D(Coordinates,NNodes,L,GN_Handle);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Add Neumann contributions onto right hand side
   
  L(NNodes(1)) = L(NNodes(1)) - FHandle(Coordinates(NNodes(1)),NNodes(1),varargin{:});
  if(size(NNodes,1) > 1)  
    L(NNodes(2)) = L(NNodes(2)) - FHandle(Coordinates(NNodes(2)),NNodes(2),varargin{:}); 
  end
  
return