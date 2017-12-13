function Lloc = LOAD_Curl_Bnd_Tangent_DGLFE2(Edge,Normal,BdFlag,Data,QuadRule,s,SHandle, ...
                               FHandle1,FHandle2,varargin)
% LOAD_CURL_BND_TANGENT_DGLFE2 Element load vector for tangential component of boundary load data.
%
%   LLOC = LOAD_CURL_BND_TANGENT_DGLFE2(EDGE,NORMAL,BDFLAG,DATA,QUADRULE,S,SHANDLE,FHANDLE)
%   computes the entries of the element load vector for the tangential component of the boundary 
%   load data.
%
%   EDGE is 2-by-2 matrix whose rows contain the start and end node of the
%   current edge.
%
%   NORMAL is 1-by-2 marix which contains the unit normal with respect to
%   the current edge EDGE.
%
%   The integer BDFLAG denotes the boundary flag of the current edge. Note
%   that for interior edges only values larger than are allowed.
%
%   The struct DATA contains the left or right hand side element data:
%    ELEMENT  Integer specifying the neighbouring element.
%    ELEMFLAG Integer specifying the element flag of the neighbouring
%             element or zero.
%    VERTICES 3-by-2 matrix specifying the vertices of the neighbouring
%             element.
%    EDGELOC  Integer specifying the local edge number on the neighbouring
%             element.
%    MATCH    Integer specifying the relative orientation of the edge with
%             respect to the orientation of the neighbouring element.
%
%   The integer S can specifies wheter the diffusive fluxes are discretized
%   in a symmetric or anti-symmetric way:
%    +1 Antisymmetric discretization of diffusive fluxes
%    -1 Symmetric discretization of diffusive fluxes
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%
%   SHANDLE is a function pointer to the edge weight function.
%
%   FHANDLE is a function pointer to the load data.
%
%   LLOC = LOAD_CURL_BND_TANGENT_DGLFE2(EDGE,NORMAL,BDFLAG,DATA,QUADRULE,SHANDLE,FHANDLE, ...
%   PARAM) also handles the variable length argumet list PARAM to the
%   function pointers SHANDLE and FHANDLE.


%   2010-2010 Chak Shing Lee
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialze constants
  
  nPts = size(QuadRule.x,1);

  % Preallocate memory
  
  Lloc = zeros(6,1);
  
  % Compute value of jump weight function
  
  P0 = Edge(1,:);
  P1 = Edge(2,:);
  sigma = SHandle(P0,P1,varargin{:});
 
  % Compute values of shape functions
  
  shap = shap_DGLFE([QuadRule.x zeros(nPts,1)]);
  
  bK = Data.Vertices(1,:);
  BK = [Data.Vertices(2,:)-bK; ...
        Data.Vertices(3,:)-bK];
  area = abs(det(BK));
  
  N = zeros(nPts,3);
  dN = zeros(1,6);
  switch(Data.EdgeLoc)
    case 1
      if(Data.Match == 1)
        N(:,1) = shap(:,3);
        N(:,2) = shap(:,1);
        N(:,3) = shap(:,2);  
      else
        N(:,1) = shap(:,3);
        N(:,2) = shap(:,2);
        N(:,3) = shap(:,1);
      end
    case 2
      if(Data.Match == 1)
        N(:,1) = shap(:,2);
        N(:,2) = shap(:,3);
        N(:,3) = shap(:,1);
      else
        N(:,1) = shap(:,1);
        N(:,2) = shap(:,3);
        N(:,3) = shap(:,2);
      end
    case 3
      if(Data.Match == 1)
        N(:,1) = shap(:,1);
        N(:,2) = shap(:,2);
        N(:,3) = shap(:,3);  
      else
        N(:,1) = shap(:,2);
        N(:,2) = shap(:,1);
        N(:,3) = shap(:,3);
      end
  end
  
  dN = -[ Data.Vertices(3,:) - Data.Vertices(2,:) ...
        Data.Vertices(1,:) - Data.Vertices(3,:) ...
        Data.Vertices(2,:) - Data.Vertices(1,:) ]/(area);
  
  % Compute element map
  
  x = ones(nPts,1)*P0 + QuadRule.x*(P1-P0);
  dS = norm(P1-P0);

  
  % Compute function values
  
  FVal1 = FHandle1(x,BdFlag,varargin{:});
  FVal2 = FHandle2(x,BdFlag,varargin{:});
  
  % Compute the tangential component 
  
  FVal = -Normal(2)*FVal1+Normal(1)*FVal2;
  
  % Compute entries of element load vector

  Lloc(1) = -Normal(2)*sigma*sum(QuadRule.w.*FVal.*N(:,1))*dS ...
            + s*dN(1)*sum(QuadRule.w.*FVal)*dS;
  Lloc(2) = Normal(1)*sigma*sum(QuadRule.w.*FVal.*N(:,1))*dS ...
            + s*dN(2)*sum(QuadRule.w.*FVal)*dS;      
  Lloc(3) = -Normal(2)*sigma*sum(QuadRule.w.*FVal.*N(:,2))*dS ...
            + s*dN(3)*sum(QuadRule.w.*FVal)*dS;
  Lloc(4) = Normal(1)*sigma*sum(QuadRule.w.*FVal.*N(:,2))*dS ...
            + s*dN(4)*sum(QuadRule.w.*FVal)*dS;      
  Lloc(5) = -Normal(2)*sigma*sum(QuadRule.w.*FVal.*N(:,3))*dS ...
            + s*dN(5)*sum(QuadRule.w.*FVal)*dS;
  Lloc(6) = Normal(1)*sigma*sum(QuadRule.w.*FVal.*N(:,3))*dS ...
            + s*dN(6)*sum(QuadRule.w.*FVal)*dS;      
  
return