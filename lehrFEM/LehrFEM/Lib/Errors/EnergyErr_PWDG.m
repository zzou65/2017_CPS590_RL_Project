function err = EnergyErr_PWDG(Mesh,u,QuadRule,omega,UHandle,DUHandle,eta,varargin)
%ENERGYERR_PWDG Helmholtz energy-norm discretization error for PWDG
% 
%   ERR = ENERGYERR_PWDG(MESH,U,QUADRULE,OMEGA,UHANDLE,DUHANDLE,ETA)
%   computes the Helmholtz energy-norm discretization error of the 
%   discontinuous plane wave finite element solution given by the vector U 
%   compared to the exact solution given by the exact solution given by the
%   function handles UHANDLE and DUHANDLE.
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
%   UHANDLE is a function handle for the exact solution and DUHANDLE is a
%   function handle for its gradient, returned row-wise.
%
%   The scalar OMEGA is the wave number of the plane waves.
%
%   ETA is the weight for the jump terms.  A value of 0 corresponds to the
%   broken Sobolev norm; a value of 1 corresponds to the standard DG norm.
%
%   ERR = EnERGYERR_PWDG(...,PARAMS) passes the variable length parameter
%   list PARAMS to the function handles UHANDLE and DUHANDLE.  These should
%   take the arguments X (spacial coordinates), element flag, PARAMS{:}.
%
%   See also L2Err_PWDG, STIMA_Lapl_Vol_PWDG, STIMA_norm_Inn_PWDG.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Intialize constants
  nPts = size(QuadRule.w,1);
  nElements = size(Mesh.Elements,1);
  nVert = size(Mesh.Elements,2);
  nDofs = [Mesh.ElemData.nDofs];    % Number of degrees of freedom on element
  nDofsSum = cumsum([0,nDofs]);     % Number of degrees of freedom on all preceding elements

  % Compute volume terms of discretization error

  err = 0;
  for j = 1:nElements
       
    % Extract vertex and basis function numbers
    vidx = Mesh.Elements(j,:);
    idx = nDofsSum(j) + (1:nDofs(j));
    
    % If all edges of element are straight / not parametrized
    if(~isfield(Mesh.ElemData,'Geom') || ~Mesh.ElemData(j).Geom.isParam)
          
      if(nVert==3) % triangles

        % Compute element mapping    
        bK = Mesh.Coordinates(vidx(1),:);
        BK = [Mesh.Coordinates(vidx(2),:)-bK; Mesh.Coordinates(vidx(3),:)-bK];
        det_BK = abs(det(BK));

        % Transform quadrature points
        x_loc = QuadRule.x*BK;
        x = x_loc + ones(nPts,1)*bK;

      elseif(nVert==4) % quadrilaterals

        % Compute element mapping
        P1 = Mesh.Coordinates(vidx(1),:);
        P2 = Mesh.Coordinates(vidx(2),:);
        P3 = Mesh.Coordinates(vidx(3),:);
        P4 = Mesh.Coordinates(vidx(4),:);

        N_BFE = shap_BFE(QuadRule.x);
        grad_N_BFE = grad_shap_BFE(QuadRule.x);

        z1 = P1(1)*grad_N_BFE(:,1:2)+P2(1)*grad_N_BFE(:,3:4) + ...
             P3(1)*grad_N_BFE(:,5:6)+P4(1)*grad_N_BFE(:,7:8);
        z2 = P1(2)*grad_N_BFE(:,1:2)+P2(2)*grad_N_BFE(:,3:4) + ...
             P3(2)*grad_N_BFE(:,5:6)+P4(2)*grad_N_BFE(:,7:8);
        det_BK = abs(z1(:,1).*z2(:,2)-z1(:,2).*z2(:,1));

        % Transform quadrature points
        x = N_BFE(:,1)*P1+N_BFE(:,2)*P2+N_BFE(:,3)*P3+N_BFE(:,4)*P4;
        x_loc = x - Mesh.Coordinates(vidx(ones(size(x,1),1)),:);

      else % not triangles or quadrilaterals
        error('This code only supports triangles and quadrilaterals');

      end

      % Evaluate basis functions
      Dir = Mesh.ElemData(j).Dir;
      N = exp(i*omega*x_loc*Dir.');

      % Evaluate solutions
      u_FE = N*u(idx);
      u_EX = UHandle(x,varargin{:});

      % Evaluate gradients of solutions
      Du_FE = i*omega*[(N.*Dir(:,ones(nPts,1)).')*u(idx),(N.*Dir(:,2*ones(nPts,1)).')*u(idx)];
      Du_EX = DUHandle(x,varargin{:});

      % Compute error on current element
      err = err + omega^2*sum(QuadRule.w.*abs(u_EX-u_FE).^2.*det_BK) ...
        + sum(QuadRule.w.*sum(abs((Du_EX-Du_FE).^2),2).*det_BK);
      
    % If an edge of the current element is paramtrized
    else
      
      % Initialize lower order quadrature rule
      QR0 = P7O6();
      QR1 = P3O3();
      
      % Load triangular mesh on element
      Mesh0 = Mesh.ElemData(j).Geom.LocMesh;
      X0 = Mesh.Coordinates(Mesh.Elements(j,1),:);
      
      % Find maximal refinement level
      maxRef = max(Mesh0.ElemFlag);
      midRef = ceil(maxRef/2);
      
      % Loop over elements of submesh and compute contributions
      for l=1:size(Mesh0.Elements,1)
        
        % Choose quadrature rule
        if(Mesh0.ElemFlag(l)>=maxRef)
          QR = QR1;
        elseif(Mesh0.ElemFlag(l)>=midRef)
          QR = QR0;
        else
          QR = QuadRule;
        end
        nQR = size(QR.x,1);

        % Get vertices of element
        Vert0 = Mesh0.Coordinates(Mesh0.Elements(l,:),:);

        % Compute element mapping
        bK = Vert0(1,:);
        BK = [Vert0(2,:) - bK; ...
              Vert0(3,:) - bK];
        det_BK = abs(det(BK));

        x_loc = QR.x*BK;
        x = x_loc + ones(nQR,1)*bK;
        x0 = bK - X0;

        % Evaluate shape functions
        Dir = Mesh.ElemData(j).Dir;
        N = exp(i*omega*(x_loc + x0(ones(nQR,1),:))*Mesh.ElemData(j).Dir.');
        
        % Evaluate solutions
        u_FE = N*u(idx);
        u_EX = UHandle(x,varargin{:});

        % Evaluate gradients of solutions
        Du_FE = i*omega*[(N.*Dir(:,ones(nQR,1)).')*u(idx),(N.*Dir(:,2*ones(nQR,1)).')*u(idx)];
        Du_EX = DUHandle(x,varargin{:});

        % Compute error on current element
        for k = 1:Mesh.ElemData(j).nDofs
          err = err + omega^2*sum(QR.w.*abs(u_EX-u_FE).^2.*det_BK) ...
            + sum(QR.w.*sum(abs((Du_EX-Du_FE).^2),2).*det_BK);
        end

      end
      
      % Clear auxiliary mesh
      clear Mesh0;
      
    end
    
  end
  
  % Handle jump terms
  if(eta)
  
    % Assemble inner product for jump terms
    B = assemMat_Inn_PDG2(Mesh,@STIMA_norm_Inn_PWDG,omega);

    % Add jump terms to error
    err = err + eta*abs(u'*B*u);
    
  end
  
  % Calculate actual error: square root of above sum
  err = sqrt(err);
  
return