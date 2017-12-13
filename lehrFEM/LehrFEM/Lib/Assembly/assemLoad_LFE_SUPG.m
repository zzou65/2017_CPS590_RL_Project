function L = assemLoad_LFE_SUPG(Mesh,QuadRule,VHandle,FHANDLE,e,d1,d2,varargin)
% ASSEMLOAD_LFE_SUPG Assemble linear FE contributions for SUPG method.
%
%   L = ASSEMLOAD_LFE_SUPG(MESH,QUADRULE,VHANDLE,FHANDLE) assembles the global load 
%   vector (SUPG-modification) for the load data given by the function handle FHANDLE.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%
%   e: diffusivity
%   d1, d2: apriori chosen constants for SUPG-modification
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.
%
%   VHANDLE is function handle for velocity field
%
%   FHANDLE is function handle for right hand side function handle
%
%   Copyright 2005-2007 Patrick Meury, Holger Heumann, Alan Mitchel
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
    P1=Mesh.Coordinates(vidx(1),:);
    P2=Mesh.Coordinates(vidx(2),:);
    P3=Mesh.Coordinates(vidx(3),:);
    
    bK = P1;
    BK = [P2-bK; P3-bK];
    inv_BK=inv(BK).';
    det_BK = abs(det(BK));
    
    N_loc(:,1:2)=N(:,1:2)*inv_BK; 
    N_loc(:,3:4)=N(:,3:4)*inv_BK;
    N_loc(:,5:6)=N(:,5:6)*inv_BK;
    
    x = QuadRule.x*BK + ones(nPts,1)*bK;
    c = VHandle(x,varargin{:});
     
    % local PecletNumber
    hK=max([norm(P2-P1),norm(P3-P1),norm(P2-P3)]);
    v_infK=max(abs(c(:)));
    PK=v_infK*hK/(2*e);
    
    % scalarproduct of velocity and Grad Lambda
    
    c1=sum(c.*[N_loc(:,1) N_loc(:,2)],2);
    c2=sum(c.*[N_loc(:,3) N_loc(:,4)],2);
    c3=sum(c.*[N_loc(:,5) N_loc(:,6)],2);
     
    % Compute load data
    dummy=0;
    FVal = FHANDLE(x,dummy,e,varargin{:});
    
    % Add contributions to global load vector
    if (PK<=1)
      d=d1*hK^2/e;
    else
      d=d2*hK;
    end
    
    L(vidx(1)) = L(vidx(1)) + d*sum(QuadRule.w.*FVal(:,1).*c1)*det_BK;
    L(vidx(2)) = L(vidx(2)) + d*sum(QuadRule.w.*FVal(:,1).*c2)*det_BK;
    L(vidx(3)) = L(vidx(3)) + d*sum(QuadRule.w.*FVal(:,1).*c3)*det_BK;
      
  end
  
return