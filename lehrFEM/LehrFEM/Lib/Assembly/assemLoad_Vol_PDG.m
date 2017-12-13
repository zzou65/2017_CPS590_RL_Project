function L = assemLoad_Vol_PDG(Mesh,nDofs,EHandle,varargin)
% ASSEMLOAD_VOL_PDG Assemble discontinuous volume contributions.
%
%   L = ASSEMLOAD_VOL_PDG(MESH,NDOFS,EHANDLE) assembles the global load
%   vector from the local element contributions given by the function
%   handle EHANDLE. The integer NDOFS specifies the number of dofs per
%   element.
%
%   A = ASSEMLOAD_VOL_PDG(MESH,NDOFS,EHANDLE,EPARAM) handles the variable
%   length argument list EPARAM to the function handle EHANDLE during the
%   assembly process. 
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%
%   Example:
%
%   F = @(x,varargin)-4*ones(size(x,1),1);
%   L = assemLoad_Vol_PDG(Mesh,3,@LOAD_Vol_PDG,P3O3(),@shap_DGCR,F);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nElements = size(Mesh.Elements,1);  % Number of elements
  
  % Preallocate memory
  
  L = zeros(nDofs*nElements,1);
  
  % Check for element flags
  
  if(isfield(Mesh,'ElemFlag')),
    ElemFlag = Mesh.ElemFlag; 
  else
    ElemFlag = zeros(nElements,1);
  end
  BdFlags = Mesh.BdFlags;
  
  % Assemble element contributions
  
  for i = 1:nElements
    
    % Extract vertices
    
    vidx = Mesh.Elements(i,:);
    
    % Compute element contributions
    
    Lloc = EHandle(Mesh.Coordinates(vidx,:),ElemFlag(i),varargin{:});
    
    % Add contributions to global load vector
    
    idx = nDofs*(i-1)+(1:nDofs);
    
    L(idx) = L(idx)+Lloc;
    
  end
  
return