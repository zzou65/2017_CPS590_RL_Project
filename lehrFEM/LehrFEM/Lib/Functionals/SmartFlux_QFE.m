function J = SmartFlux_QFE(u,Mesh,PsiHandle,QuadRule,SigmaHandle,varargin)
% SMARTFLUX_QFE Flux throu specified part of the boundary.
%
%   J = SMARTFLUX_QFE(U,MESH,PSIHANDLE,QUADRULE,SIGMAHANDLE) computes the
%   value of the flux throu the boundary specified by the cut-off function
%   PSIHANDLE of the finite element solution U on the struct MESH using
%   the quadrature rule QUADRULE.
% 
%   J = SMARTFLUX_QFE(U,MESH,PSIHANDLE,QUADRULE,SIGMAHANDLE,SIGMAPARAM)
%   also handles the variable length argument list SIGMAPARAM to the
%   function handle SIGMAPARAM during assembly.
%
%   Example:
%
%   PsiHandle = @(x,varargin)[x(:,1) x(:,2)];
%   SigmaHandle = @(x,varargin)ones(size(x,1),1);
%   J = SmartFlux_QFE(u,Mesh,PsiHandle,QuadRule,SigmaHandle)

%   Copyright 2005-2005 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nElements = size(Mesh.Elements,1);
  nCoordinates = size(Mesh.Coordinates,1);
  nGauss = size(QuadRule.w,1);
  
  % Precompute gradients of shape functions
  
  grad_N = grad_shap_QFE(QuadRule.x);
  
  J = 0;
  for i = 1:nElements
     
    % Extract vertices of current element
    
    vidx = Mesh.Elements(i,:);
    Vertices = Mesh.Coordinates(vidx,:);
    
    % Compute element mapping
    
    bK = Vertices(1,:);
    BK = [Vertices(2,:)-bK; Vertices(3,:)-bK];
    inv_BK_t = transpose(inv(BK));
    det_BK = abs(det(BK));
    
    x = QuadRule.x*BK+ones(nGauss,1)*bK;
 
    % Compute function values on current element
    
    SigmaVal = SigmaHandle(x,Mesh.ElemFlag(i),varargin{:});
    grad_Psi = PsiHandle(x,Mesh.ElemFlag(i)); 
    idx = [vidx ...
           Mesh.Vert2Edge(vidx(1),vidx(2))+nCoordinates ...
           Mesh.Vert2Edge(vidx(2),vidx(3))+nCoordinates ...
           Mesh.Vert2Edge(vidx(3),vidx(1))+nCoordinates];
    grad_U = (u(idx(1))*grad_N(:,1:2)+u(idx(2))*grad_N(:,3:4) + ...
              u(idx(3))*grad_N(:,5:6)+u(idx(4))*grad_N(:,7:8) + ...
              u(idx(5))*grad_N(:,9:10)+u(idx(6))*grad_N(:,11:12))*inv_BK_t;
    
    % Update flux 
    
    J = J + sum(QuadRule.w.*SigmaVal.*sum(grad_U.*grad_Psi,2))*det_BK;
      
  end
  
return