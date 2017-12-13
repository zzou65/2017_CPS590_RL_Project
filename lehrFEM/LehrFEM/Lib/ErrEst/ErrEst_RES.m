function Eta = ErrEst_RES(U,Mesh,QuadRule,FHandle,varargin)
% ERREST_RES Residual based error estimator for the Laplacian.
%
%   ETA = ERREST_RES(U,MESH,QUADRULE,FHANDLE) computes the residual based
%   error estimator ETA of the solution U to the Dirichlet problem for the 
%   Laplace equation.
%   
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    EDGE2ELEM    P-by-2 matrix connecting edges to elements. The first 
%                 column specifies the left hand side element where the
%                 second column specifies the right hand side element.
%
%   Example:
%
%   Eta = ErrEst_RES(U,Mesh,P7O6(),gauleg(0,1,3),F_HANDLE);
 
%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nElements = size(Mesh.Elements,1);
  nEdges = size(Mesh.Edges,1);
  nGauss = size(QuadRule.w,1);
  
  % Preallocate memory
  
  Eta = zeros(nElements,1);
   
  % Compute element contributions
  
  for i = 1:nElements
          
    vidx = Mesh.Elements(i,:);
    P1 = Mesh.Coordinates(vidx(1),:);
    P2 = Mesh.Coordinates(vidx(2),:);
    P3 = Mesh.Coordinates(vidx(3),:);
    
    % Compute the element width
    
    h_K = max([norm(P1-P2) norm(P2-P3) norm(P3-P1)]);
    
    % Compute element mapping
    
    bK = P1;
    BK = [ P2 - P1 ; P3 - P1 ]; 
    inv_BK = inv(BK);
    det_BK = abs(det(BK));
    x = QuadRule.x*BK + ones(nGauss,1)*bK;
    
    % Compute error contribution on each element
    
    FVal = FHandle(x,varargin{:});      
    Eta(i) = h_K^2*sum(QuadRule.w.*abs(FVal).^2)*det_BK;
    
  end

  % Compute edge contributions
  
  grad_N = grad_shap_LFE([0 0]);
  Rot = [0 -1; 1 0];
  for i = 1:nEdges
    if(Mesh.BdFlags(i) >= 0)        
      P0 = Mesh.Coordinates(Mesh.Edges(i,1),:);
      P1 = Mesh.Coordinates(Mesh.Edges(i,2),:);
        
      % Compute unit normal and edge length
      
      normal = P1-P0;
      h_F = sqrt(sum(normal.^2));
      normal = normal*Rot/h_F;
        
      % Compute left and right hand side neighbours  
        
      Elem_l = Mesh.Edge2Elem(i,1);
      vidx_l = Mesh.Elements(Elem_l,:);
      Elem_r = Mesh.Edge2Elem(i,2);
      vidx_r = Mesh.Elements(Elem_r,:);
      
      % Compute element mappings
      
      bK_l = Mesh.Coordinates(vidx_l(1),:);
      BK_l = [Mesh.Coordinates(vidx_l(2),:)-bK_l; ...
              Mesh.Coordinates(vidx_l(3),:)-bK_l];
      bK_r = Mesh.Coordinates(vidx_r(1),:);
      BK_r = [Mesh.Coordinates(vidx_r(2),:)-bK_r; ...
              Mesh.Coordinates(vidx_r(3),:)-bK_r];
          
      inv_BK_l = inv(BK_l);
      inv_BK_r = inv(BK_r);
            
      % Compute left and right hand-side gradients
      
      grad_u_l = (U(vidx_l(1))*grad_N(1:2) + ...
                  U(vidx_l(2))*grad_N(3:4) + ...
                  U(vidx_l(3))*grad_N(5:6))*transpose(inv_BK_l);
      grad_u_r = (U(vidx_r(1))*grad_N(1:2) + ...
                  U(vidx_r(2))*grad_N(3:4) + ...
                  U(vidx_r(3))*grad_N(5:6))*transpose(inv_BK_r);
      
      % Add edge error contributions to left and right hand-side neighbours
      
      Eta_F = h_F^2*abs(sum((grad_u_l-grad_u_r).*normal,2))^2;   
      Eta(Elem_l) = Eta(Elem_l) + Eta_F/2;
      Eta(Elem_r) = Eta(Elem_r) + Eta_F/2;
      
    end
  end
  
  Eta = sqrt(Eta);
  
return