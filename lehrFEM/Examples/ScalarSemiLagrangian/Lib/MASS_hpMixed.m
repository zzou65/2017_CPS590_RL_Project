function Mloc = MASS_hpMixed(ElemInfo,EDofs,EDir,CDofs,Shap, p_EDofs,p_EDir,p_CDofs,p_Shap,weights,varargin)
% MASS_HP Element mass matrix.
%
%   MLOC = MASS_HPMixed(ElemInfo,EDofs,EDir,CDofs,Shap, p_EDofs,p_EDir,p_CDofs,p_Shap,weights,
%   computes the value of the element mass matrix using hierarchical shape
%   functions.
%
%   ELEMINFO is an integer parameter which is used to specify additional
%   element information on each element.
%
%   EDOFS and P_EDOFS are a 1-by-3 matrix specifying the maximum polynomial
%   degree of
%   the edge shape functions on every edge.
%
%   EDIR and P_EDIR are 1-by-3 matrix specifying the orientation of the edges of the
%   element.
%
%   CDOFS  and P_CDOFS are an integer specifying the maximum polynomial degree inside the
%   element.
%
%   QUADRULE is a struct, which specifies a 2D Gauss qaudrature that is
%   used to do the integration on each element:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.
%
%   The struct SHAP and p_SHAP contains the values and gradients of the shape
%   functions computed by the routine SHAP_HP at the quadrature points of
%   QUADRULE.
%   
%   WEIGHTS contains the scaled quadrature weights of QUADRULE
%
%   Example:
%
%   p = 10;
%   q =10;
%   EDofs = (p-1)*ones(1,3);
%   EDir = [1 1 1];
%   CDofs = (p-1)*(p-2)/2;
%   p_EDofs = (q-1)*ones(1,3);
%   p_EDir = [1 1 1];
%   p_CDofs = (q-1)*(q-2)/2;
%   QuadRule = Duffy(TProd(gauleg(0,1,2*p)));
%   Shap = shap_hp(QuadRule.x,p);
%   p_Shap = shap_hp(QuadRule.x,q);
%   weights = QuadRule.w;
%   Mloc2 = MASS_hpMixed(0,EDofs,EDir,CDofs,Shap,p_EDofs,p_EDir,p_CDofs,p_Shap, weights);

%   Copyright 2006-2011 Patrick Meury and Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Preallocate memory
  
  Mloc = zeros(3+sum(EDofs)+CDofs,3+sum(p_EDofs)+p_CDofs);

  for i = 1:3
      shap_1 = Shap.vshap.shap{i};

      % Vertex vs vertex
      for j = 1:3
          shap_2 = p_Shap.vshap.shap{j};
          Mloc(i,j) = sum( weights.*shap_1.*shap_2) ;
      end

      % Vertex vs edge 1
      offset = 3;
      for j = 1:p_EDofs(1)
          shap_2 = p_EDir(1)^(j+1)*p_Shap.eshap{1}.shap{j};
          Mloc(i,j+offset) = sum( weights.*shap_1.*shap_2);
      end

      % Vertex vs edge 2
      offset = 3+p_EDofs(1);
      for j = 1:p_EDofs(2)
          shap_2 = p_EDir(2)^(j+1)*p_Shap.eshap{2}.shap{j};
          Mloc(i,j+offset) = sum( weights.*shap_1.*shap_2) ;
      end

      % Vertex vs edge 3
      offset = 3+p_EDofs(1)+p_EDofs(2);
      for j = 1:p_EDofs(3)
          shap_2 = p_EDir(3)^(j+1)*p_Shap.eshap{3}.shap{j};
          Mloc(i,j+offset) = sum( weights.*shap_1.*shap_2) ;
      end

      % Vertex vs element
      offset = 3+sum(p_EDofs);
      for j = 1:p_CDofs
          shap_2 = p_Shap.cshap.shap{j};
          Mloc(i,j+offset) = sum( weights.*shap_1.*shap_2) ;
      end
  end

  offset_1 = 3;
  for i = 1:EDofs(1)
      shap_1 = EDir(1)^(i+1)*Shap.eshap{1}.shap{i};
      
      % Edge1 vs vertex
      for j = 1:3
          shap_2 = p_Shap.vshap.shap{j};
          Mloc(i+offset_1,j) = sum( weights.*shap_1.*shap_2) ;
      end

      % Edge 1 vs edge 1
      offset_2 = 3;
      for j = 1:p_EDofs(1)
          shap_2 = p_EDir(1)^(j+1)*p_Shap.eshap{1}.shap{j};
          Mloc(i+offset_1,j+offset_2) = sum( weights.*shap_1.*shap_2) ;
      end

      % Edge 1 vs edge 2
      offset_2 = 3+p_EDofs(1);
      for j = 1:p_EDofs(2)
          shap_2 = p_EDir(2)^(j+1)*p_Shap.eshap{2}.shap{j};
          Mloc(i+offset_1,j+offset_2) = sum( weights.*shap_1.*shap_2) ;
      end

      % Edge 1 vs edge 3
      offset_2 = 3+p_EDofs(1)+p_EDofs(2);
      for j = 1:p_EDofs(3)
          shap_2 = p_EDir(3)^(j+1)*p_Shap.eshap{3}.shap{j};
          Mloc(i+offset_1,j+offset_2) = sum( weights.*shap_1.*shap_2) ;
      end

      % Edge 1 vs element
      offset_2 = 3+sum(p_EDofs);
      for j = 1:p_CDofs
          shap_2 = p_Shap.cshap.shap{j};
          Mloc(i+offset_1,j+offset_2) = sum( weights.*shap_1.*shap_2) ;
      end
  end

  offset_1 = 3+EDofs(1);
  for i = 1:EDofs(2)
      shap_1 = EDir(2)^(i+1)*Shap.eshap{2}.shap{i};
      
      % Edge2 vs vertex
      for j = 1:3
          shap_2 = p_Shap.vshap.shap{j};
          Mloc(i+offset_1,j) = sum( weights.*shap_1.*shap_2) ;
      end

      % Edge 2 vs edge 1
      offset_2 = 3;
      for j = 1:p_EDofs(1)
          shap_2 = p_EDir(1)^(j+1)*p_Shap.eshap{1}.shap{j};
          Mloc(i+offset_1,j+offset_2) = sum( weights.*shap_1.*shap_2) ;
      end

      % Edge 2 vs edge 2
      offset_2 = 3+p_EDofs(1);
      for j = 1:p_EDofs(2)
          shap_2 = p_EDir(2)^(j+1)*p_Shap.eshap{2}.shap{j};
          Mloc(i+offset_1,j+offset_2) = sum( weights.*shap_1.*shap_2) ;
      end

      % Edge 2 vs edge 3
      offset_2 = 3+p_EDofs(1)+p_EDofs(2);
      for j = 1:p_EDofs(3)
          shap_2 = p_EDir(3)^(j+1)*p_Shap.eshap{3}.shap{j};
          Mloc(i+offset_1,j+offset_2) = sum( weights.*shap_1.*shap_2) ;
      end

      % Edge 2 vs element
      offset_2 = 3+sum(p_EDofs);
      for j = 1:p_CDofs
          shap_2 = p_Shap.cshap.shap{j};
          Mloc(i+offset_1,j+offset_2) = sum( weights.*shap_1.*shap_2) ;
      end
  end

  offset_1 = 3+EDofs(1)+EDofs(2);
  for i = 1:EDofs(3)
      shap_1 = EDir(3)^(i+1)*Shap.eshap{3}.shap{i};

      % Edge 3 vs vertex
      for j = 1:3
          shap_2 = p_Shap.vshap.shap{j};
          Mloc(i+offset_1,j) = sum( weights.*shap_1.*shap_2) ;
      end

      % Edge 3 vs edge 1
      offset_2 = 3;
      for j = 1:p_EDofs(1)
          shap_2 = p_EDir(1)^(j+1)*p_Shap.eshap{1}.shap{j};
          Mloc(i+offset_1,j+offset_2) = sum( weights.*shap_1.*shap_2) ;
      end

      % Edge 3 vs edge 2
      offset_2 = 3+p_EDofs(1);
      for j = 1:p_EDofs(2)
          shap_2 = p_EDir(2)^(j+1)*p_Shap.eshap{2}.shap{j};
          Mloc(i+offset_1,j+offset_2) = sum( weights.*shap_1.*shap_2) ;
      end
      
      % Edge 3 vs edge 3
      offset_2 = 3+p_EDofs(1)+p_EDofs(2);
      for j = 1:p_EDofs(3)
          shap_2 = p_EDir(3)^(j+1)*p_Shap.eshap{3}.shap{j};
          Mloc(i+offset_1,j+offset_2) = sum( weights.*shap_1.*shap_2) ;
      end

      % Edge 3 vs element
      offset_2 = 3+sum(p_EDofs);
      for j = 1:p_CDofs
          shap_2 = p_Shap.cshap.shap{j};
          Mloc(i+offset_1,j+offset_2) = sum( weights.*shap_1.*shap_2) ;
      end
  end

  offset_1 = 3+sum(EDofs);
  for i = 1:CDofs
      shap_1 = Shap.cshap.shap{i};

      % Element vs vertex
      for j = 1:3
          shap_2 = p_Shap.vshap.shap{j};
          Mloc(i+offset_1,j) = sum( weights.*shap_1.*shap_2) ;
      end

      % Element vs edge 1
      offset_2 = 3;
      for j = 1:p_EDofs(1)
          shap_2 = p_EDir(1)^(j+1)*p_Shap.eshap{1}.shap{j};
          Mloc(i+offset_1,j+offset_2) = sum( weights.*shap_1.*shap_2) ;
      end

      % Element vs edge 2
      offset_2 = 3+p_EDofs(1);
      for j = 1:p_EDofs(2)
          shap_2 = p_EDir(2)^(j+1)*p_Shap.eshap{2}.shap{j};
          Mloc(i+offset_1,j+offset_2) = sum( weights.*shap_1.*shap_2) ;
      end
      
      % Element vs edge 3
      offset_2 = 3+p_EDofs(1)+p_EDofs(2);
      for j = 1:p_EDofs(3)
          shap_2 = p_EDir(3)^(j+1)*p_Shap.eshap{3}.shap{j};
          Mloc(i+offset_1,j+offset_2) = sum( weights.*shap_1.*shap_2) ;
      end
      
      % Element vs element
      offset_2 = 3+sum(p_EDofs);
      for j = 1:p_CDofs
          shap_2 = p_Shap.cshap.shap{j};
          Mloc(i+offset_1,j+offset_2) = sum( weights.*shap_1.*shap_2) ;
      end
  end
return