function Eta = ErrEst_REC(U,Mesh,QuadRule)
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
  nGauss = size(QuadRule.w,1);
  
  % Compute recovered gradient
    
  grad_N = grad_shap_LFE([0 0]);
  grad_U_REC = zeros(nCoordinates,2);  
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
      
      Gu = Gu + (U(vidx(1))*grad_N(1:2) + ...
                 U(vidx(2))*grad_N(3:4) + ...
                 U(vidx(3))*grad_N(5:6))*transpose(inv_BK)*det_BK/2;                   
      area = area + det_BK/2;
      
    end
    grad_U_REC(i,:) = Gu/area;
  end
  
  % Compute error estimator
  
  N = shap_LFE(QuadRule.x);
  grad_N = grad_shap_LFE(QuadRule.x);
  Eta = zeros(nElements,1); 
  for i = 1:nElements
    vidx = Mesh.Elements(i,:);
      
    % Compute element mapping
    
    bK = Mesh.Coordinates(vidx(1),:);
    BK = [Mesh.Coordinates(vidx(2),:)-bK; ...
          Mesh.Coordinates(vidx(3),:)-bK];
    inv_BK = inv(BK);
    det_BK = abs(det(BK));
    
    % Evaluate gradient and recovered gradient at quadrature points
    
    grad_u = (U(vidx(1))*grad_N(:,1:2) + ...
              U(vidx(2))*grad_N(:,3:4) + ...
              U(vidx(3))*grad_N(:,5:6))*transpose(inv_BK);
    Gu = N(:,1)*grad_U_REC(vidx(1),:) + ...
         N(:,2)*grad_U_REC(vidx(2),:) + ...
         N(:,3)*grad_U_REC(vidx(3),:);
  
    % Compute error estimate on the current element
     
    Eta(i) = sqrt(sum(QuadRule.w.*sum((grad_u-Gu).^2,2))*det_BK);
    
  end
  
return
