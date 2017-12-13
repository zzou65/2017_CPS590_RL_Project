function Aloc = STIMA_Lapl_hp(Vertices,ElemInfo,EDofs,EDir,CDofs,QuadRule,Shap,varargin)
% STIMA_LAPL_HP Element stiffness matrix for the Laplacian.
%
%   ALOC = STMA_LAPL_HP(VERTICES,ELEMINFO,EDOFS,EDIR,CDOFS,QUADRULE,SHAP)
%   computes the value of the element stiffness matrix for the Laplacian
%   using hierarchical shape functions.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current
%   element in a row wise orientation.
%
%   ELEMINFO is an integer parameter which is used to specify additional
%   element information on each element.
%
%   EDOFS is a 1-by-3 matrix specifying the maximum polynomial degree of
%   the edge shape functions on every edge.
%
%   EDIR is 1-by-3 matrix specifying the orientation of the edges of the
%   element.
%
%   CDOFS is an integer specifying the maximum polynomial degree inside the
%   element.
%
%   QUADRULE is a struct, which specifies a 2D Gauss qaudrature that is
%   used to do the integration on each element:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.
%
%   The struct SHAP contains the values and gradients of the shape
%   functions computed by the routine SHAP_HP at the quadrature points of
%   QUADRULE.
%
%   Example:
%
%   p = 10;
%   Vertices = [0 0; 1 0; 0 1];
%   EDofs = (p-1)*ones(1,3);
%   EDir = [1 1 1];
%   CDofs = (p-1)*(p-2)/2;
%   QuadRule = Duffy(TProd(gauleg(0,1,2*p)));
%   Shap = shap_hp(QuadRule.x,p);
%   Aloc = STIMA_Lapl_hp(Vertices,0,EDofs,EDir,CDofs,QuadRule,Shap);

%   Copyright 2006-2006 Patrick Meury
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
   % Preallocate memory
  
  Aloc = zeros(3+sum(EDofs)+CDofs,3+sum(EDofs)+CDofs);
  
  % Compute element map
  
  bK = Vertices(1,:);
  BK = [Vertices(2,:)-bK; Vertices(3,:)-bK];
  inv_BK = inv(BK);
  det_BK = abs(det(BK));
  TK = det_BK*transpose(inv_BK)*inv_BK;
  
  for i = 1:3
    grad_shap_1 = Shap.vshap.grad_shap{i};  
      
    % Vertex vs vertex
      
    for j = i:3
      grad_shap_2 = Shap.vshap.grad_shap{j};
      Aloc(i,j) = sum(QuadRule.w.*sum(grad_shap_1.*(grad_shap_2*TK),2));
    end
      
    % Vertex vs edge 1
    
    offset = 3;
    for j = 1:EDofs(1)
      grad_shap_2 = EDir(1)^(j+1)*Shap.eshap{1}.grad_shap{j};
      Aloc(i,j+offset) = sum(QuadRule.w.*sum(grad_shap_1.*(grad_shap_2*TK),2));
    end
    
    % Vertex vs edge 2
    
    offset = 3+EDofs(1);
    for j = 1:EDofs(2)
      grad_shap_2 = EDir(2)^(j+1)*Shap.eshap{2}.grad_shap{j};
      Aloc(i,j+offset) = sum(QuadRule.w.*sum(grad_shap_1.*(grad_shap_2*TK),2));    
    end
    
    % Vertex vs edge 3
    
    offset = 3+EDofs(1)+EDofs(2);
    for j = 1:EDofs(3)
      grad_shap_2 = EDir(3)^(j+1)*Shap.eshap{3}.grad_shap{j};
      Aloc(i,j+offset) = sum(QuadRule.w.*sum(grad_shap_1.*(grad_shap_2*TK),2));
    end
    
    % Vertex vs element
    
    offset = 3+sum(EDofs);
    for j = 1:CDofs
      grad_shap_2 = Shap.cshap.grad_shap{j};
      Aloc(i,j+offset) = sum(QuadRule.w.*sum(grad_shap_1.*(grad_shap_2*TK),2));     
    end
    
  end
           
  offset_1 = 3;
  for i = 1:EDofs(1)
    grad_shap_1 = EDir(1)^(i+1)*Shap.eshap{1}.grad_shap{i};  
      
    % Edge 1 vs edge 1  
      
    offset_2 = 3;
    for j = i:EDofs(1)
      grad_shap_2 = EDir(1)^(j+1)*Shap.eshap{1}.grad_shap{j};  
      Aloc(i+offset_1,j+offset_2) = sum(QuadRule.w.*sum(grad_shap_1.*(grad_shap_2*TK),2));  
    end
    
    % Edge 1 vs edge 2
    
    offset_2 = 3+EDofs(1);
    for j = 1:EDofs(2)
      grad_shap_2 = EDir(2)^(j+1)*Shap.eshap{2}.grad_shap{j};
      Aloc(i+offset_1,j+offset_2) = sum(QuadRule.w.*sum(grad_shap_1.*(grad_shap_2*TK),2));
    end
    
    % Edge 1 vs edge 3
    
    offset_2 = 3+EDofs(1)+EDofs(2);
    for j = 1:EDofs(3)       
      grad_shap_2 = EDir(3)^(j+1)*Shap.eshap{3}.grad_shap{j};
      Aloc(i+offset_1,j+offset_2) = sum(QuadRule.w.*sum(grad_shap_1.*(grad_shap_2*TK),2));   
    end
      
    % Edge 1 vs element
    
    offset_2 = 3+sum(EDofs);
    for j = 1:CDofs
      grad_shap_2 = Shap.cshap.grad_shap{j};
      Aloc(i+offset_1,j+offset_2) = sum(QuadRule.w.*sum(grad_shap_1.*(grad_shap_2*TK),2));  
    end
    
  end
  
  offset_1 = 3+EDofs(1);
  for i = 1:EDofs(2)
    grad_shap_1 = EDir(2)^(i+1)*Shap.eshap{2}.grad_shap{i};  
      
    % Edge 2 vs edge 2  
      
    offset_2 = 3+EDofs(1);
    for j = i:EDofs(2)
      grad_shap_2 = EDir(2)^(j+1)*Shap.eshap{2}.grad_shap{j};
      Aloc(i+offset_1,j+offset_2) = sum(QuadRule.w.*sum(grad_shap_1.*(grad_shap_2*TK),2));  
    end
    
    % Edge 2 vs edge 3
    
    offset_2 = 3+EDofs(1)+EDofs(2);
    for j = 1:EDofs(3)
      grad_shap_2 = EDir(3)^(j+1)*Shap.eshap{3}.grad_shap{j};
      Aloc(i+offset_1,j+offset_2) = sum(QuadRule.w.*sum(grad_shap_1.*(grad_shap_2*TK),2));
    end    
       
    % Edge 2 vs element
    
    offset_2 = 3+sum(EDofs);
    for j = 1:CDofs
      grad_shap_2 = Shap.cshap.grad_shap{j};
      Aloc(i+offset_1,j+offset_2) = sum(QuadRule.w.*sum(grad_shap_1.*(grad_shap_2*TK),2));
    end    
        
  end
  
  offset_1 = 3+EDofs(1)+EDofs(2);
  for i = 1:EDofs(3)
    grad_shap_1 = EDir(3)^(i+1)*Shap.eshap{3}.grad_shap{i};  
      
    % Edge 3 vs edge 3
      
    offset_2 = 3+EDofs(1)+EDofs(2);
    for j = i:EDofs(3)
      grad_shap_2 = EDir(3)^(j+1)*Shap.eshap{3}.grad_shap{j};
      Aloc(i+offset_1,j+offset_2) = sum(QuadRule.w.*sum(grad_shap_1.*(grad_shap_2*TK),2));  
    end
     
    % Edge 3 vs element
    
    offset_2 = 3+sum(EDofs);
    for j = 1:CDofs
      grad_shap_2 = Shap.cshap.grad_shap{j};
      Aloc(i+offset_1,j+offset_2) = sum(QuadRule.w.*sum(grad_shap_1.*(grad_shap_2*TK),2));  
    end
    
  end
  
  % Element vs element
  
  offset = 3+sum(EDofs);
  for i = 1:CDofs
    grad_shap_1 = Shap.cshap.grad_shap{i};      
    for j = i:CDofs
      grad_shap_2 = Shap.cshap.grad_shap{j};
      Aloc(i+offset,j+offset) = sum(QuadRule.w.*sum(grad_shap_1.*(grad_shap_2*TK),2));     
    end
  end
   
  % Fill in lower trinagular part
    
  tri = triu(Aloc);
  Aloc = tril(transpose(tri),-1)+tri;
  
return
