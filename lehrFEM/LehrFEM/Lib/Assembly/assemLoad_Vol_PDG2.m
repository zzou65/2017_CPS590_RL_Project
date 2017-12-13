function L = assemLoad_Vol_PDG2(Mesh,EHandle,varargin)
%ASSEMLOAD_VOL_PDG2 Assemble DG volume contributions
%
%   L = ASSEM_LOAD_VOL_PDG2(MESH,EHANDLE) assembles the global load vector
%   from the local element contributions given by the function handle
%   EHANDLE.
%
%   L = ASSEM_LOAD_VOL_PDG2(MESH,EHANDLE,EPARAM) passes the variable-length
%   argument list EPARAM to the function handle EHANDLE during the assembly
%   process.
%
%   The struct MESH must contain at least the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    ELEMDATA     N-by-1 structure array containing at least the fields:
%       NDOFS       The number of degrees of freedom on the corresponding
%                   element.
%
%   See also assemLoad_Bnd_PDG2, assemMat_Vol_PDG2, assemMat_Inn_PDG2,
%   assemMat_Bnd_PDG2.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  nElements = size(Mesh.Elements,1); % Number of elements

  % Count degrees of freedom on all elements and calculate number of
  % matrix entries
  nDofs = [Mesh.ElemData.nDofs];    % Number of degrees of freedom on element
  nDofsSum = cumsum([0,nDofs]);     % Number of degrees of freedom on all preceding elements
  numel = nDofsSum(end);            % Number of vector entries
  
  % Preallocate memory
  L = zeros(numel,1);
  
  % Assemble element contributions
  for i = 1:nElements
    
    % Extract vertices
    vidx = Mesh.Elements(i,:);
    
    % Compute element contributions
    Lloc = EHandle(Mesh.Coordinates(vidx,:),Mesh.ElemData(i),varargin{:});
    
    % Add contributions to global load vector
    idx = nDofsSum(i)+(1:nDofs(i));
    L(idx) = L(idx)+Lloc;
    
  end
  
return