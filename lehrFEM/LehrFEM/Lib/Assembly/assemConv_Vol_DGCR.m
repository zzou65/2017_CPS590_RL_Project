function B = assemConv_Vol_DGCR(Mesh,U,Lim,QuadRule,FHandle,varargin)
% ASSEMCONV_VOL_DGCR Assemble non-linear volume convection terms.
%
%   B = ASSEMCON_VOL_DGCR(MESH,U,LIM,QUADRULE,FHANDLE) assembles the global
%   load vector from the local element contributions.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%
%   U denotes the FE ansatz function for which the volume convection term
%   has to be assembled.
%
%   The integer LIM denotes the type of limiter to be used for the test
%   functions:
%    0  Piecwise constant test functions 
%    1  Piecewise linear test functions
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%   
%   The function handle FHANDLE denotes the non-linear flux function use
%   for the computation of the volume convection term.
%    
%   B = ASSEMCONV_VOL_DGCR(MESH,U,LIM,QUADRULE,FHANDLE,PARAM) also handles
%   the variable length argument list PARAM to the flux function FHANDLE.
%
%   Example:
%
%   B = assemConv_Vol_DGCR(Mesh,U,1,P3O3(),@f_LA);
%
%   See also grad_shap_DGCR.

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nElements = size(Mesh.Elements,1);
  nPts = size(QuadRule.x,1);
  
  % Preallocate memory
  
  B = zeros(3*nElements,1);
  
  % Check for limited test functions
  
  if(Lim > 0)
  
    % Check for element flags
  
    if(isfield(Mesh,'ElemFlag')),
      ElemFlag = Mesh.ElemFlag; 
    else
      ElemFlag = zeros(nElements,1);
    end
  
    % Precompute values and gradients of shape functions
  
    N = shap_DGCR(QuadRule.x);
    grad_shap = grad_shap_DGCR(QuadRule.x);
  
    % Assemble convection contributions
  
    for i = 1:nElements
        
      % Extract vertex numbers
    
      vidx = Mesh.Elements(i,:);
    
      % Compute element mapping
    
      bK = Mesh.Coordinates(vidx(1),:);
      BK = [Mesh.Coordinates(vidx(2),:)-bK; ...
            Mesh.Coordinates(vidx(3),:)-bK];    
      inv_BK_t = transpose(inv(BK));
      det_BK = abs(det(BK));
    
      x = QuadRule.x*BK + ones(nPts,1)*bK;
           
      % Compute function value of non-linear part 
    
      idx = 3*(i-1)+[1 2 3];
      UVal = U(idx(1))*N(:,1)+U(idx(2))*N(:,2)+U(idx(3))*N(:,3);
      FVal = FHandle(x,UVal,ElemFlag(i),varargin{:});

      % Add convection contributions 
          
      grad_N = grad_shap(:,1:2)*inv_BK_t;
      B(idx(1)) = sum(QuadRule.w.*sum(FVal.*grad_N,2))*det_BK;
    
      grad_N = grad_shap(:,3:4)*inv_BK_t;
      B(idx(2)) = sum(QuadRule.w.*sum(FVal.*grad_N,2))*det_BK;
    
      grad_N = grad_shap(:,5:6)*inv_BK_t;
      B(idx(3)) = sum(QuadRule.w.*sum(FVal.*grad_N,2))*det_BK;
     
    end
    
  end

return
