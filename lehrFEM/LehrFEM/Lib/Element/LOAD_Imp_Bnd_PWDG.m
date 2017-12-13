function Lloc = LOAD_Imp_Bnd_PWDG(Edge,Normal,Params,Data,QuadRule,omega,GHandle,varargin)
%LOAD_IMP_BND_PWDG Element load vector for impedance boundary conditions
%
%   LLOC = LOAD_IMP_BND_PWDG(EDGE,NORMAL,PARAMS,DATA,QUADRULE,OMEGA,GHANDLE)
%   computes the enrtries of the element load vector corresponding to 
%   impedance boundary conditions for discontinuous plane waves.
%
%   EDGE is 2-by-2 matrix whose rows contain the start and end node of the
%   current edge.
%
%   NORMAL is 1-by-2 marix which contains the interior unit normal with
%   respect to the current edge EDGE.
%
%   PARAMS is a structure that contains at least the fields:
%    D        Scalar coefficient for flux on impedence boundary, should be
%             between 0 and 1.
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
%   GHANDLE is a function handle for the impedence boundary conditions.  It
%   takes as arguments the spacial coordinates and the exterior normal
%   vector.
%
%   LLOC = LOAD_IMP_BND_PWDG(...,GPARAM) passes the variable length
%   parameter list GPARAM to the function handle GHANDLE.
%
%   Example:
%
%     B_Vol = assemLoad_Vol_PDG2(Mesh,@LOAD_Vol_PWDG,qr2,omega,f);
%     B_Imp = assemLoad_Bnd_PDG2(Mesh,-1,@LOAD_Imp_Bnd_PWDG,qr1,omega,gI);
%     B_Dir = assemLoad_Bnd_PDG2(Mesh,-2,@LOAD_Dir_Bnd_PWDG,qr1,omega,gD);
%     B = B_Vol + B_Imp + B_Dir;
%
%   See also assemLoad_Bnd_PDG2, LOAD_Dir_Bnd_PWDG, LOAD_Neu_Bnd_PWDG,
%   LOAD_Vol_PWDG, STIMA_Helm_Imp_Bnd_PWDG, STIMA_Helm_Inn_PWDG.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize
  nDofs = Data.ElemData.nDofs;      % Number of degrees of freedom
  Dir = Data.ElemData.Dir;          % Propagation directions of plane wave basis functions
  
  % Check for parametrization of edge
  isParam = isfield(Params,'Geom') && Params.Geom.isParam;
  
  % If edge is straight / not parametrized
  if(~isParam)

    % Orient normal vector as exterior normal vector
    n = -Normal;
  %   orient = check_ext_normal(Normal,Edge,Data.Vertices,Data.EdgeLoc);
  %   n = orient*Normal;

    % The fluxes will be represented as linear combinations of plane waves
    % with direction vectors Dir and value one at Edge(1,:)

    % Calculate traces of test vectors (normal components)
    d = Edge(1,:) - Data.Vertices(1,:);
    v = exp(i*omega*Dir*d.');
    Dv = i*omega*v.*(Dir*n.');

    % Calculate plane wave test function components
    v0 = (1-Params.d)*conj(v) - Params.d/(i*omega)*conj(Dv);

    % Calculate remaining integral
    ba = Edge(2,:) - Edge(1,:);
    nx = numel(QuadRule.x);
    Dirba = Dir*ba.';
    Lloc = v0*sqrt(sum(ba.^2));
    GVal = GHandle(Edge(ones(nx,1),:)+ba(ones(nx,1),:).*QuadRule.x(:,[1 1]),n,varargin{:});
    for j=1:nDofs
      Lloc(j) = Lloc(j)*sum(QuadRule.w.*exp(-i*omega*Dirba(j)*QuadRule.x).*GVal);
    end
  
  % If edge is parametrized
  else
    
    % Calculate integral directly through quadrature
    nx = numel(QuadRule.x);
    X0 = Data.Vertices(1,:);
    Gamma = Params.Geom.Gamma(QuadRule.x);
    N = Params.Geom.N(QuadRule.x);
    GVal = GHandle(Gamma,N,varargin{:});
    Lloc = Params.Geom.dGamma(ones(nDofs,1));
    for j=1:nDofs
      ExpVal = exp(i*omega*(Gamma-X0(ones(nx,1),:))*Dir(j,:).');
      VVal = (1-Params.d)*conj(ExpVal) - (Params.d/(i*omega))*conj(i*omega*ExpVal.*(N*Dir(j,:).'));
      Lloc(j) = Lloc(j)*sum(QuadRule.w.*VVal.*GVal);
    end
    
  end
  
return