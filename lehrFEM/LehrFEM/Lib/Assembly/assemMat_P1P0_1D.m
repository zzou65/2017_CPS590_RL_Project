function A = assemMat_P1P0_1D(Coordinates,EHandle,varargin)
% ASSEMMAT_P1P0_1D Assemble linear/constant FE contributions.
%
%   A = ASSEMMAT_P1P0_1D(COORDINATES,EHANDLE) assembles the global matrix
%   from the local element contributions given by the function handle
%   EHANDLE and returns the matrix in a sparse representation.
%
%   A = ASSEMMAT_P1P0_1D(COORDINATES,EHANDLE,EPARAM) handles the variable
%   length argument list EPARAM to the function handle EHANDLE during the
%   assembly process. 
%
%   Example:
%
%   Coordinates = transpose(0:.01:1);
%   A = assemMat_P1P0_1D(Coordinates,@STIMA_Lapl_1D);

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nCoordinates = size(Coordinates,1);
  nElements = nCoordinates-1;
  
  % Preallocate memory
  
  Diag = zeros(nCoordinates,2);
  
  vidx = [1 2];
  for i = 1:nElements

    % Compute element stiffness matrix
      
    Aloc = EHandle(Coordinates(vidx),varargin{:});  
    
    % Add contributions to global matrix
    
    Diag(vidx(1),2) = Diag(vidx(1),2) + Aloc(1);
    Diag(vidx(1),1) = Diag(vidx(1),1) + Aloc(2);
        
    % Update current element
    
    vidx = vidx+1;
    
  end
  
  % Assign output arguments
  
  A = spdiags(Diag,-1:0,nCoordinates,nElements);    
  
return