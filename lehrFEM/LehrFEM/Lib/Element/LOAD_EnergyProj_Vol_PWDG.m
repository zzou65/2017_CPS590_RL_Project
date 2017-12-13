function Lloc = LOAD_EnergyProj_Vol_PWDG(Vertices,ElemData,QuadRule,omega,UHandle,DUHandle,varargin)
%LOAD_ENERGYPROJ_VOL_PWDG element load vector for energy-norm projection
%
%   LLOC = LOAD_ENERGYPROJ_VOL_PWDG(VERTICES,ELEMDATA,QUADRULE,OMEGA,
%   UHANDLE,DUHANDLE) computes the volume contributions to the element load
%   vector for discontinuous plane waves with right hand side given by an
%   energy-norm scalar product.
%
%   VERTICES is 3-by-2 or 4-by-2 matrix whose rows contain the vertices of
%   the current element.
%
%   ELEMDATA is a structure that contains at least the fields:
%     NDOFS   The number of degrees of freedom on the current element.
%     DIR     A P-by-2 matrix containing the propagation directions of
%             the plane wave basis functions in its rows.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%
%   The scalar OMEGA is the wave number of the plane waves.
%
%   UHANDLE and DUHANDLE are function handles for the load data and its
%   gradient; they should take at least the argument x.  To pass further
%   arguments to these function handles, append them to the end of the
%   list of input arguments.
%
%   See also assemLoad_Vol_PDG2, LOAD_Vol_PWDG, EnergyErr_PWDG.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialze constants
  nPts = size(QuadRule.x,1);  % Number of quadrature points
  nVert = size(Vertices,1);   % Number of edges of polygon
  nDofs = ElemData.nDofs;     % Number of basis functions
  Dir = ElemData.Dir;         % Propagation directions of plane wave basis functions
  
  % Check for parametrization of edge
  isParam = isfield(ElemData,'Geom') && ElemData.Geom.isParam;
  
  % If all edges are straight / not parametrized
  if(~isParam)
  
    % Compute element mapping
    if(nVert==3) % triangle
      bK = Vertices(1,:);
      BK = [Vertices(2,:) - bK; ...
            Vertices(3,:) - bK];
      det_BK = abs(det(BK));

      x_loc = QuadRule.x*BK;
      x = x_loc + ones(nPts,1)*bK;

    elseif(nVert==4) % quadrilateral   
      P1 = Vertices(1,:);
      P2 = Vertices(2,:);
      P3 = Vertices(3,:);
      P4 = Vertices(4,:);

      N_BFE = shap_BFE(QuadRule.x);
      grad_N_BFE = grad_shap_BFE(QuadRule.x);

      z1 = P1(1)*grad_N_BFE(:,1:2)+P2(1)*grad_N_BFE(:,3:4) + ...
           P3(1)*grad_N_BFE(:,5:6)+P4(1)*grad_N_BFE(:,7:8);
      z2 = P1(2)*grad_N_BFE(:,1:2)+P2(2)*grad_N_BFE(:,3:4) + ...
           P3(2)*grad_N_BFE(:,5:6)+P4(2)*grad_N_BFE(:,7:8);
      det_BK = abs(z1(:,1).*z2(:,2)-z1(:,2).*z2(:,1));

      x = N_BFE(:,1)*P1+N_BFE(:,2)*P2+N_BFE(:,3)*P3+N_BFE(:,4)*P4;
      x_loc = x - Vertices(ones(size(x,1),1),:);

    else % not triangle or quadrilateral
      error('This code only supports triangles and quadrilaterals');

    end

    % Evaluate shape functions
    N = exp(i*omega*x_loc*Dir.');

    % Compute function values
    UVal = UHandle(x,varargin{:});
    DUVal = DUHandle(x,varargin{:});

    % Preallocate memory
    Lloc = zeros(nDofs,1);

    % Compute entries of element load vector
    for j = 1:nDofs
      Lloc(j) = omega^2*sum(QuadRule.w.*UVal.*conj(N(:,j)).*det_BK) ...
        + sum(QuadRule.w.*sum(DUVal.*conj(i*omega*Dir(j*ones(nPts,1),:).*N(:,[j j])),2).*det_BK);
    end
    
  % If an edge of the element is parametrized
  else
    
    % Load triangular mesh on element
    Mesh = ElemData.Geom.LocMesh;
      
    % Initialize lower order quadrature rule
    QR0 = P7O6();
    QR1 = P3O3();
      
    % Find maximal refinement level
    maxRef = max(Mesh.ElemFlag);
    midRef = ceil(maxRef/2);
    
    % Initialize output
    Lloc = zeros(ElemData.nDofs,1);
    X0 = Vertices(1,:);
    
    % Loop over elements of fine mesh and compute contributions
    for j=1:size(Mesh.Elements,1)
      
      % Choose quadrature rule
      if(Mesh.ElemFlag(j)>=maxRef)
        QR = QR1;
      elseif(Mesh.ElemFlag(j)>=midRef)
        QR = QR0;
      else
        QR = QuadRule;
      end
      nQR = size(QR.x,1);
      
      % Get vertices of element
      Vert0 = Mesh.Coordinates(Mesh.Elements(j,:),:);
      
      % Compute element mapping
      bK = Vert0(1,:);
      BK = [Vert0(2,:) - bK; ...
            Vert0(3,:) - bK];
      det_BK = abs(det(BK));

      x_loc = QR.x*BK;
      x = x_loc + ones(nQR,1)*bK;
      x0 = bK - X0;

      % Evaluate shape functions
      N = exp(i*omega*(x_loc + x0(ones(nQR,1),:))*ElemData.Dir.');

      % Compute function values
      UVal = UHandle(x,varargin{:});
      DUVal = DUHandle(x,varargin{:});

      % Compute entries of element load vector
      for k = 1:ElemData.nDofs
        Lloc(k) = Lloc(k) + omega^2*sum(QR.w.*UVal.*conj(N(:,k)).*det_BK) ...
          + sum(QR.w.*sum(DUVal.*conj(i*omega*Dir(k*ones(nQR,1),:).*N(:,[k k])),2).*det_BK);
      end
      
    end
    
  end
  
return