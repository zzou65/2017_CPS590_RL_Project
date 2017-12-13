function Lloc = LOAD_Dir_Bnd_PWDG_dual(Edge,Normal,Params,Data,QuadRule,omega,GHandle,varargin)
%LOAD_DIR_BND_PWDG_DUAL Element load vector for Dirichlet boundary cond.
%
%   LLOC =LOAD_DIR_BND_PWDG(EDGE,NORMAL,PARAMS,DATA,QUADRULE,OMEGA,GHANDLE)
%   computes the entries of the element load vector corresponding to 
%   Dirichlet boundary conditions for dual (mixed) discontinuous plane
%   waves.
%
%   EDGE is 2-by-2 matrix whose rows contain the start and end node of the
%   current edge.
%
%   NORMAL is 1-by-2 marix which contains the interior unit normal with
%   respect to the current edge EDGE.
%
%   PARAMS is a structure that contains at least the fields:
%    A        Scalar coefficient for the (jump-jump) penalty term in the
%             numerical flux for the function value.
%    NDOFS    Total number of degrees of freedom on elements adjacent to
%             current edge.
%
%   The struct DATA contains the left or right hand side element data:
%    ELEMENT  Integer specifying the neighbouring element.
%    ELEMDATA Structure containing at least the fields:
%       NDOFS   The number of degrees of freedom on the current element.
%       DIR     A P-by-2 matrix containing the propagation directions of
%               the plane wave basis functions in its rows.
%    VERTICES 3-by-2 or 4-by-2 matrix specifying the vertices of the
%             neighbouring element.
%    EDGELOC  Integer specifying the local edge number on the neighbouring
%             element.
%    MATCH    Integer specifying the relative orientation of the edge with
%             respect to the orientation of the neighbouring element.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    w Weights of the Gauss quadrature.
%    x Abscissae of the Gauss quadrature.
%   
%   OMEGA is the wave number of the Helholtz equation.
%
%   GHANDLE is a function handle for the Dirichlet boundary conditions.  It
%   takes as arguments the spacial coordinates.
%
%   LLOC = LOAD_IMP_BND_PWDG_DUAL(...,GPARAM) passes the variable length
%   parameter list GPARAM to the function handle GHANDLE.
%
%   Example:
%
%     B_Vol = assemLoad_Vol_PDG2_vec(...
%       Mesh,3,@LOAD_Vol_PWDG_dual,qr2,omega,f);
%     B_Imp = assemLoad_Bnd_PDG2_vec(...
%       Mesh,-1,3,@LOAD_Imp_Bnd_PWDG_dual,qr1,omega,gI);
%     B_Dir = assemLoad_Bnd_PDG2_vec(...
%       Mesh,-2,3,@LOAD_Dir_Bnd_PWDG_dual,qr1,omega,gD);
%     B = B_Vol + B_Imp + B_Dir;

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize
  nDofs = Data.ElemData.nDofs;      % Number of degrees of freedom
  Dir = Data.ElemData.Dir;          % Propagation directions of plane wave basis functions

  % Orient normal vector as exterior normal vector
  n = -Normal;
%   orient = check_ext_normal(Normal,Edge,Data.Vertices,Data.EdgeLoc);
%   n = orient*Normal;
  
  % The fluxes will be represented as linear combinations of plane wave
  % basis functions with direction vectors Dir.  Note that, since
  % the values are three-dimensional vectors (function & gradient), this
  % representation will have dimension 3*nDofs.

  % Constuct normal vector in vector form
  N = [zeros(nDofs,1);n(ones(nDofs,1),1);n(ones(nDofs,1),2)];
  
  % Calculate traces (normal components)
  d = Edge(1,:) - Data.Vertices(1,:);
  v0 = exp(i*omega*Dir*d.');
  v = [v0;zeros(2*nDofs,1)];
  tau = [zeros(nDofs,1);v0;v0].*N;
  
  % Calculate plane wave test function components
  v0 = Params.a*conj(v) + conj(tau);
  
  % Calculate remaining integral
  ba = Edge(2,:) - Edge(1,:);
  nx = numel(QuadRule.x);
  Dirba = Dir*ba.';
  Lloc = v0*sqrt(sum(ba.^2));
  GVal = GHandle(Edge(ones(nx,1),:)+ba(ones(nx,1),:).*QuadRule.x(:,[1 1]),varargin{:});
  ind = (0:2)*nDofs;
  for j=1:nDofs
    Lloc(j+ind) = Lloc(j+ind)*sum(QuadRule.w.*exp(-i*omega*Dirba(j)*QuadRule.x).*GVal);
  end
  
return