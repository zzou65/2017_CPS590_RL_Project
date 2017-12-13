function L = assemLoad_PML_Vol(Mesh,QuadRule,EpsHandle,FHANDLE,varargin)
% ASSEMLOAD_PML_Vol.
%
%   L = ASSEMLOAD_PML_Vol(MESH,QUADRULE,VHANDLE,FHANDLE) assembles the global load 
%   vector for the load data given by the function handle FHANDLE.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.
%
%   EpsHANDLE is function handle for tenor weight
%
%   FHANDLE is function handle for right hand side function handle
%
%   L = ASSEMLOAD_LFE(MESH,QUADRULE,VHandle,FHANDLE,FPARAM) also handles the 
%   additional variable length argument list FPARAM to the function handle
%   FHANDLE.
%
%   Example:
%
%   FHandle = @(x,varargin)x(:,1).^2+x(:,2).^2;
%   VHandle = @(x,varargin)[x(:,1) -x(:,2)];
%   L = assemLoad_PML_Vol(Mesh,P7O6(),epsHandle,FHandle);
%
%   See also grad_shap_LFE.

%   Copyright 2005-2007 Patrick Meury, Holger Heumann, Alan Mitchel, Nana
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(QuadRule.w,1);
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  
  % Preallocate memory
  
  L = zeros(nCoordinates,1);
  
  % Precompute shape functions
  
  N = grad_shap_LFE(QuadRule.x);
   
  % Assemble element contributions
  
  for i = 1:nElements
    
    % Extract vertices
    
    vidx = Mesh.Elements(i,:);  
    
    % Compute element mapping
    
    bK = Mesh.Coordinates(vidx(1),:);
    BK = [Mesh.Coordinates(vidx(2),:)-bK; Mesh.Coordinates(vidx(3),:)-bK];
    inv_BK=inv(BK).';
    det_BK = abs(det(BK));
    
    N_loc(:,1:2)=N(:,1:2)*inv_BK; 
    N_loc(:,3:4)=N(:,3:4)*inv_BK;
    N_loc(:,5:6)=N(:,5:6)*inv_BK;
    
    x = QuadRule.x*BK + ones(nPts,1)*bK;
    EVal = epsHandle(x);
    
    FVal = FHANDLE(x);
    
    % scalarproduct of velocity and Grad Lambda
    
    Val=eVal(:,1)*Fval(:,1)*N_loc(:,1)+...
        eVal(:,2)*Fval(:,2)*N_loc(:,1)+...
        eVal(:,3)*Fval(:,1)*N_loc(:,2)+...
        eVal(:,4)*Fval(:,2)*N_loc(:,2);
     
    % Compute load data
    
    % Add contributions to global load vector
    
    L(vidx(1)) = L(vidx(1)) + sum(QuadRule.w.*FVal(:,1).*c1)*det_BK;
    L(vidx(2)) = L(vidx(2)) + sum(QuadRule.w.*FVal(:,1).*c2)*det_BK;
    L(vidx(3)) = L(vidx(3)) + sum(QuadRule.w.*FVal(:,1).*c3)*det_BK;
      
  end
  
return