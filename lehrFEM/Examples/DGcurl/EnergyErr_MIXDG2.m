function err = EnergyErr_MIXDG2(Mesh,u,QuadRule,UHandle1,UHandle2,Curl_UHandle,SIGMA1,SIGMA2,alpha1,alpha2,varargin)
%EnergyErr_MIXDG2 DG energy-norm of discretization error of the discontinuous galerkin solution
% 
%   ERR = ENERGYERR_MIXDG2(Mesh,u,QuadRule,UHANDLE1,UHANDLE2,Curl_UHANDLE,SIGMA1,SIGMA2,alpha1,alpha2,varargin)
%   computes the DG energy-norm of discretization error of the 
%   discontinuous galerkin solution given by the vector U 
%   compared to the exact solution given by the
%   function handles UHANDLE1 and UHANDLE2.
%
%   The struct MESH should contain at least the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 or N-by-4 matrix specifying the elements of the
%                mesh.
%    ELEMDATA    N-by-1 structure array containing at least the fields:
%       NDOFS    	 The number of degrees of freedom on the corresponding
%                  element.
%       DIR        A NDOFS-by-2 matrix containing the propagation
%                  directions of the plane wave basis functions in
%                  its rows.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.
%
%   UHANDLE1 is a function handle for the first component of the exact solution and UHANDLE2 is a
%   function handle for the second component of the exact solution.
%
%   Curl_UHANDLE is function handles for the curl of the exact solution.
%
%   SIGMA1 and SIGMA2 are function handles for the stabilization term in
%   the normal and tangential jump, respectively.

%   2010-2010 Chak Shing Lee
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants
  nPts = size(QuadRule.w,1);
  nElements = size(Mesh.Elements,1);
  nCoordinates = size(Mesh.Coordinates,1);
  nEdges = size(Mesh.Edges,1);
  %nVert = size(Mesh.Elements,2);
  %nDofs = [Mesh.ElemData.nDofs];    % Number of degrees of freedom on element
  %nDofsSum = cumsum([0,nDofs]);     % Number of degrees of freedom on all preceding elements

  % Compute volume terms of discretization error

  
  if(isfield(Mesh,'ElemFlag')),
    ElemFlag = Mesh.ElemFlag; 
  else
    ElemFlag = zeros(nElements,1);
  end
  % Precompute shape function values at the quedrature points
  
  N = shap_DGLFE(QuadRule.x);
  
  err = 0;
  for i = 1:nElements
       
    % Extract vertex and basis function numbers
    vidx = Mesh.Elements(i,:);
    idx = 6*(i-1)+[1 2 3 4 5 6];

    bK = Mesh.Coordinates(vidx(1),:);
    BK = [Mesh.Coordinates(vidx(2),:)-bK; Mesh.Coordinates(vidx(3),:)-bK];
    det_BK = abs(det(BK));
    
    dNC = -[ Mesh.Coordinates(vidx(3),:) - Mesh.Coordinates(vidx(2),:) ...
        Mesh.Coordinates(vidx(1),:) - Mesh.Coordinates(vidx(3),:) ...
        Mesh.Coordinates(vidx(2),:) - Mesh.Coordinates(vidx(1),:) ]/det_BK;
    
    % Transform quadrature points
      
    x = QuadRule.x*BK+ones(nPts,1)*bK;
      
    % Evaluate solutions
      
    u_EX1 = UHandle1(x,ElemFlag(i),varargin{:});
    u_EX2 = UHandle2(x,ElemFlag(i),varargin{:});
    u_FE1 = u(idx(1))*N(:,1) + u(idx(3))*N(:,2) + u(idx(5))*N(:,3);
    u_FE2 = u(idx(2))*N(:,1) + u(idx(4))*N(:,2) + u(idx(6))*N(:,3);

    % Evaluate curl and divergence of solutions
    Curl_u_FE = u(idx(1))*dNC(1) + u(idx(2))*dNC(2) + u(idx(3))*dNC(3) + u(idx(4))*dNC(4) + u(idx(5))*dNC(5) + u(idx(6))*dNC(6);
    Curl_u_EX = Curl_UHandle(x,varargin{:});

    % Compute error on current element
    err = err+sum(QuadRule.w.*abs(u_EX1-u_FE1).^2)*det_BK+sum(QuadRule.w.*abs(u_EX2-u_FE2).^2)*det_BK+sum(QuadRule.w.*abs(Curl_u_EX-Curl_u_FE).^2)*det_BK;
    
  end
  
 % Handle interior jump terms
 [I2,J2,Jinn_n] = assemMat_Inn_DG2(Mesh,@STIMA_InnPen_Normal_DGLFE2,SIGMA1);
 [I2,J2,Jinn_t] = assemMat_Inn_DG2(Mesh,@STIMA_InnPen_Tangent_DGLFE2,SIGMA2);
 [I3,J3,Jbnd] = assemMat_Bnd_DG2(Mesh,@STIMA_BndPen_Tangent_DGLFE2,SIGMA2); 
 X = sparse([J2;J3], ...
            [I2;I3], ...
            [1/alpha2*Jinn_t+1/alpha1*Jinn_n;0*Jbnd]);    
 err = err+u'*X*u; 

 % Handle boundary terms
 for j = 1:nEdges
      
    if(Mesh.BdFlags(j) < 0)
        if(Mesh.Edge2Elem(j,1) > 0)
            Data.Element = Mesh.Edge2Elem(j,1);
            Data.ElemFlag = ElemFlag(Data.Element);
            Data.Vertices = Mesh.Coordinates(Mesh.Elements(Data.Element,:),:);
            Data.EdgeLoc = Mesh.EdgeLoc(j,1);
            Data.Match = Mesh.EdgeOrient(j,1);
        else
            Data.Element = Mesh.Edge2Elem(j,2);
            Data.ElemFlag = ElemFlag(Data.Element);
            Data.Vertices = Mesh.Coordinates(Mesh.Elements(Data.Element,:),:);
            Data.EdgeLoc = Mesh.EdgeLoc(j,2);
            Data.Match = Mesh.EdgeOrient(j,2);
        end
        
    Normal = Mesh.Normals(j,:);    
    QuadRule_Edge = gauleg(0,1,15);
    
    P0 = Mesh.Coordinates(Mesh.Edges(j,1),:);
    P1 = Mesh.Coordinates(Mesh.Edges(j,2),:);
    
    bK = Data.Vertices;
    BK = [bK(2,:)-bK(1,:); ...
          bK(3,:)-bK(1,:)];
    
    idx = 6*(Data.Element-1)+[1 2 3 4 5 6];
    nPts_Edge = size(QuadRule_Edge.w,1);
    x = ones(nPts_Edge,1)*P0 + QuadRule_Edge.x*(P1-P0);
    x_Ref = (x-ones(nPts_Edge,1)*bK(1,:))*inv(BK);
    N = shap_DGLFE(x_Ref);
    u_FE1 = u(idx(1))*N(:,1) + u(idx(3))*N(:,2) + u(idx(5))*N(:,3);
    u_FE2 = u(idx(2))*N(:,1) + u(idx(4))*N(:,2) + u(idx(6))*N(:,3);
    u_EX1 = UHandle1(x,ElemFlag(i),varargin{:});
    u_EX2 = UHandle2(x,ElemFlag(i),varargin{:});
    
    
    % Add jump terms to error
 %   err = err + eta*abs(u'*B*u);
    err = err + sum(QuadRule_Edge.w.*abs(Normal(2)*(u_EX1-u_FE1)-Normal(1)*(u_EX2-u_FE2)).^2);
  
    end
  end
  
  % Calculate actual error: square root of above sum
  err = sqrt(err);
  
return