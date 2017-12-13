function L = assemLoad_QFE_SUPG(Mesh,QuadRule,VHandle,FHANDLE,e,d1,d2,varargin)
% ASSEMLOAD_QFE_SUPG Assemble quadratic FE contributions for SUPG method.
%
%   L = ASSEMLOAD_QFE_SUPG(MESH,QUADRULE,VHANDLE,FHANDLE,e,d1,d2) assembles the global load 
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
%   Copyright 2005-2008 Patrick Meury, Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nPts = size(QuadRule.w,1);
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  nEdges = size(Mesh.Edges,1);
  
  % Preallocate memory
  
  L = zeros(nCoordinates+nEdges,1);
 
  % gradientents and laplacians of shape functions
  dN = grad_shap_QFE(QuadRule.x);
  ddN = lapl_shap_QFE(QuadRule.x);
   
  % Assemble element contributions
  
  for i = 1:nElements
    
    % Extract vertices
    
    vidx = Mesh.Elements(i,:);   
    eidx(1) = Mesh.Vert2Edge(vidx(1),vidx(2)) + nCoordinates;
    eidx(2) = Mesh.Vert2Edge(vidx(2),vidx(3)) + nCoordinates;
    eidx(3) = Mesh.Vert2Edge(vidx(3),vidx(1)) + nCoordinates;
    
    % Compute element mapping
    P1=Mesh.Coordinates(vidx(1),:);
    P2=Mesh.Coordinates(vidx(2),:);
    P3=Mesh.Coordinates(vidx(3),:);
    
    bK = P1;
    BK = [P2-bK; P3-bK];
    inv_BK=inv(BK);
    inv_BK_t=transpose(inv_BK);
    det_BK = abs(det(BK));

    % velocity at quadrature points
    x = QuadRule.x*BK + ones(nPts,1)*bK;
    c = VHandle(x,varargin{:});

    % gradients of shape functions
    grad_lambda=dN;
    for k=1:2:12
        grad_lambda(:,k:k+1)=grad_lambda(:,k:k+1)*inv_BK_t;
    end

    % scalarproduct of velocity and Grad Lambda
    cgl=zeros(nPts,6);
    for k=1:2:12
        cgl(:,(k+1)/2)=sum(c.*grad_lambda(:,k:k+1),2);
    end

    % Lapl of Shapfunctions
    ll_m=zeros(1,6);
    gl = grad_shap_LFE([1/3,1/3]); % constants
    gl(:,1:2)=gl(:,1:2)*inv_BK_t;
    gl(:,3:4)=gl(:,3:4)*inv_BK_t;
    gl(:,5:6)=gl(:,5:6)*inv_BK_t;
    ll_m(1)=4*gl(:,1:2)*gl(:,1:2)';
    ll_m(2)=4*gl(:,3:4)*gl(:,3:4)';
    ll_m(3)=4*gl(:,5:6)*gl(:,5:6)';
    ll_m(4)=8*gl(:,1:2)*gl(:,3:4)';
    ll_m(5)=8*gl(:,3:4)*gl(:,5:6)';
    ll_m(6)=8*gl(:,5:6)*gl(:,1:2)';
    
    % Compute load data
    dummy=0;
    FVal = FHANDLE(x,dummy,e,varargin{:});
    
    % local PecletNumber
    hK=max([norm(P2-P1),norm(P3-P1),norm(P2-P3)]);
    v_infK=max(abs(c(:)));
    PK=v_infK*hK/(2*e);
    if (PK<=1)
      d=d1*hK^2/e;
    else
      d=d2*hK;
    end
    
    % Add contributions to global load vector
    w=det_BK*QuadRule.w;
    
    L(vidx(1)) = L(vidx(1)) + d*det_BK*(...
        sum(w.*FVal.*cgl(:,1))-...
        e*ll_m(1)*sum(w.*FVal));
    L(vidx(2)) = L(vidx(2)) + d*det_BK*(...
        sum(w.*FVal.*cgl(:,2))-...
        e*ll_m(2)*sum(w.*FVal));
    L(vidx(3)) = L(vidx(3)) + d*det_BK*(...
        sum(w.*FVal.*cgl(:,3))-...
        e*ll_m(3)*sum(w.*FVal));
    L(eidx(1)) = L(eidx(1)) + d*det_BK*(...
        sum(w.*FVal.*cgl(:,4))-...
        e*ll_m(4)*sum(w.*FVal));
    L(eidx(2)) = L(eidx(2)) + d*det_BK*(...
        sum(w.*FVal.*cgl(:,5))-...
        e*ll_m(5)*sum(w.*FVal));
    L(eidx(3)) = L(eidx(3)) + d*det_BK*(...
        sum(w.*FVal.*cgl(:,6))-...
        e*ll_m(6)*sum(w.*FVal));
  end
  
return
