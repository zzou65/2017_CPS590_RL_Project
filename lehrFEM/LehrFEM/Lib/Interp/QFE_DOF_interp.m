function U = QFE_DOF_interp(Mesh,F_Handle,varargin)
% QFE_DOF_interp linear interpolation using QFE DOFs
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
%   U = QFE_DOF_interp(Mesh,F_Handle);
%
%   See also LFE_DOF_interp.

%   Copyright 2006 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  nCoordinates = size(Mesh.Coordinates,1);
  nEdges =size(Mesh.Edges,1); 
  U=zeros((nEdges+nCoordinates),1);

% Compute discretization error

  U(1:nCoordinates)=F_Handle(Mesh.Coordinates);
  eidx = zeros(1,3);
  for i = 1:nEdges
     vid=Mesh.Edges(i,:);
     U(nCoordinates+i)=F_Handle(0.5*(Mesh.Coordinates(vid(1),:)+Mesh.Coordinates(vid(2),:)));   
  end
return    
    
