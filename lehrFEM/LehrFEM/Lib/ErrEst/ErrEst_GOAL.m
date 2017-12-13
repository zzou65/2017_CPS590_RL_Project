function Eta = ErrEst_GOAL(gF,u,FHandle,Mesh,QuadRule2D,QuadRule_1D)
% ERREST_REC Recovery based error estimator.
%
%   ETA = ERREST_REC(U,MESH,QUADRULE) computes the recovery based error
%   estimator ETA of the solution U.
%   
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    ADJELEMENTS  M-by-Q matrix specifying the elements of the mesh sharing
%                 vertex i of the mesh.
%    NADJELEMENTS M-by-1 matrix specifying the actual number of elements
%                 sharing vertex COORDINATES(i) of the mesh.  
%
%   Example:
%
%   Eta = ErrEst_REC(U,Mesh,P7O6());

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  nPts = size(QuadRule2D.w,1);
  nPts_1D = size(QuadRule_1D.w,1);
  Eta = zeros(nElements,1);
  
  % Compute recovered gradient
    
  grad_N = grad_shap_LFE([0 0]);
  grad_g_REC = zeros(nCoordinates,2);  
  for i = 1:nCoordinates 
    area = 0;
    Gu = zeros(1,2);
    for j = Mesh.nAdjElements(i)
      vidx = Mesh.Elements(Mesh.AdjElements(i,j),:);
        
      % Compute element mapping
            
      bK = Mesh.Coordinates(vidx(1),:);
      BK = [Mesh.Coordinates(vidx(2),:)-bK; ...
            Mesh.Coordinates(vidx(3),:)-bK];
      inv_BK = inv(BK);
      det_BK = abs(det(BK));
      
            
      % Sum up mean value of gradient on the current patch
      
      Gu = Gu + (gF(vidx(1))*grad_N(1:2) + ...
                 gF(vidx(2))*grad_N(3:4) + ...
                 gF(vidx(3))*grad_N(5:6))*transpose(inv_BK)*det_BK/2;                   
      area = area + det_BK/2;
      
    end
    grad_g_REC(i,:) = Gu/area;
  end
  
  % Compute the H1 seminorm of the recovered gradient and error estimator
  
  grad_N = grad_shap_LFE(QuadRule2D.x);
  grad_N_1D_1 = grad_shap_LFE([QuadRule_1D.x zeros(nPts_1D,1)]);
  grad_N_1D_2 = grad_shap_LFE([zeros(nPts_1D,1) QuadRule_1D.x]);
  grad_N_1D_3 = grad_shap_LFE([QuadRule_1D.x 1-QuadRule_1D.x]);
  for i=1:nElements
    vidx = Mesh.Elements(i,:);
    a1 = Mesh.Coordinates(vidx(1),:);
    a2 = Mesh.Coordinates(vidx(2),:);
    a3 = Mesh.Coordinates(vidx(3),:);
    Gu_loc_1 = grad_g_REC(vidx,1);
    Gu_loc_2 = grad_g_REC(vidx,2);
    
    % Compute element mapping
    
    bK = Mesh.Coordinates(vidx(1),:);
    BK = [Mesh.Coordinates(vidx(2),:)-bK; ...
          Mesh.Coordinates(vidx(3),:)-bK];
    inv_BK = inv(BK);
    det_BK = abs(det(BK));
    
    grad_Gu_loc_1 = (Gu_loc_1(1)*grad_N(:,1:2)+ ...
                    Gu_loc_1(2)*grad_N(:,3:4)+ ...
                    Gu_loc_1(3)*grad_N(:,5:6))*transpose(inv_BK);
                
    grad_Gu_loc_2 = (Gu_loc_2(1)*grad_N(:,1:2)+ ...
                    Gu_loc_2(2)*grad_N(:,3:4)+ ...
                    Gu_loc_2(3)*grad_N(:,5:6))*transpose(inv_BK);
                
    H2_semi_g = sum(QuadRule2D.w.*sum(abs(grad_Gu_loc_1).^2,2))*det_BK +...
             sum(QuadRule2D.w.*sum(abs(grad_Gu_loc_2).^2,2))*det_BK;
         
    % compute quadrature points and compute norm of source term
      
    x = QuadRule2D.x*BK+ones(nPts,1)*bK;
    cell_res = sqrt(sum(QuadRule2D.w.*(FHandle(x).^2)));
    
    % compute the the L2 norm on the boundary of the current element
    
    % first edge
    
    edge = a2-a1;
    n = [edge(2)*ones(nPts_1D,1) -edge(1)*ones(nPts_1D,1)];
    grad_u = (u(vidx(1))*grad_N_1D_1(:,1:2)+ ...
             u(vidx(2))*grad_N_1D_1(:,3:4)+ ...
             u(vidx(3))*grad_N_1D_1(:,5:6))*transpose(inv_BK);
    edge1 = sum(QuadRule_1D.w .*(sum(n.*grad_u,2).^2));
    
    % second edge
    
    edge = a3-a1;
    n = [edge(2)*ones(nPts_1D,1) -edge(1)*ones(nPts_1D,1)];
    grad_u = (u(vidx(1))*grad_N_1D_2(:,1:2)+ ...
             u(vidx(2))*grad_N_1D_2(:,3:4)+ ...
             u(vidx(3))*grad_N_1D_2(:,5:6))*transpose(inv_BK);
    edge2 = sum(QuadRule_1D.w .*(sum(n.*grad_u,2).^2));
    
    % third edge
    
    edge = a3-a2;
    n = [edge(2)*ones(nPts_1D,1) -edge(1)*ones(nPts_1D,1)];
    grad_u = (u(vidx(1))*grad_N_1D_3(:,1:2)+ ...
             u(vidx(2))*grad_N_1D_3(:,3:4)+ ...
             u(vidx(3))*grad_N_1D_3(:,5:6))*transpose(inv_BK);
    edge3 = sum(QuadRule_1D.w .*(sum(n.*grad_u,2).^2));
    
    flux_res = sqrt(edge1 + edge2 + edge3);
    
    % compute the local residual
    
    h_K = max([norm(a3-a2) norm(a3-a1) norm(a1-a2)]);   
    rho_K = cell_res + 1/sqrt(h_K) * flux_res;
    
    Eta(i) = rho_K * h_K^2 * H2_semi_g;
    
  end
 
return