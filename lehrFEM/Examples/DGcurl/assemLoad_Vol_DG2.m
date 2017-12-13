function L = assemLoad_Vol_DG2(Mesh,EHandle,varargin)
% ASSEMLOAD_VOL_DG2 Assemble vectorial discontinuous Lagrangian nodal FE volume
% contributions.
%
%   L = ASSEMLOAD_VOL_DG2(MESH,EHANDLE) assembles the global load vector
%   from the local element contributions given by the function handle
%   EHANDLE.
%
%   A = ASSEMLOAD_VOL_DG2(MESH,EHANDLE,EPARAM) handles the variable length
%   argument list EPARAM to the function handle EHANDLE during the assembly
%   process. 
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.

%   2010-2010 Chak Shing Lee
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nElements = size(Mesh.Elements,1);  % Number of elements
  
  % Preallocate memory
  
  L = zeros(6*nElements,1);
  
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
    
    idx = 6*(i-1)+[1 2 3 4 5 6];
    
    L(idx) = L(idx)+Lloc;
    
  end
  
return