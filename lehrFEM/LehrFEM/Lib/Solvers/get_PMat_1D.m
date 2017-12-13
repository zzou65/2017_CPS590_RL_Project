function P = get_PMat_1D(CMesh,FMesh)
% GET_PMAT_1D Compute prolongation matrix.
%
%   P = GET_PMAT_1D(CMESH,FMESH) computes the matrix corresponding to the
%   prolongation operator from a coarse mesh CMESH to a fine mesh FMESH.
%
%   Example:
%
%   P = get_PMat_1D(CMesh,FMesh);
%
%   See also add_MLevel.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constants

    NDofs = length(CMesh)
    
% Preallocate memory

    I = zeros(3*NDofs-3,1);
    J = zeros(3*NDofs-3,1);
    M = zeros(3*NDofs-3,1);
    Contributor = zeros(3*NDofs-3,1);    

% Build prolongation matrix

    I = repmat([1 2 2]',NDofs-1,1);
    J = repmat([1 1 2]',NDofs-1,1);
    M = repmat([1 .5 .5]',NDofs-1,1);
    Contributor = repmat(0:NDofs-2,3,1);
    Contributor = Contributor(:);
    
    % Add contributor to I,J
    
    I = I + 2*Contributor;
    J = J + Contributor;
    
% Convert to sparse matrix
    
    P = sparse([I;2*NDofs-1],[J;NDofs],[M;1]);
    
 return
   