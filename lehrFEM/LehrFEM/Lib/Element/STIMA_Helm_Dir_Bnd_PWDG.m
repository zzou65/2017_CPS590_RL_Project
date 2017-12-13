function Aloc = STIMA_Helm_Dir_Bnd_PWDG(Edge,Normal,Params,Data,omega,QuadRule,varargin)
%STIMA_HELM_DIR_BND_PWDG Element stiffness matrix for PWDG Helmholtz
%
%   ALOC = STIMA_HELM_DIR_BND_PWDG(EDGE,NOTMAL,PARAMS,DATA,OMEGA)
%   computes the entries of the element stiffness matrix for a
%   discontinuous planw wave discretization of the Helmholtz equation on
%   Dirichlet boundary edges.
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
%    L2       The L2 inner product matrix on the current edge.
%
%   The structs DATA contains the left or right hand side element data:
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
%   OMEGA is the wave number of the Helholtz equation.
%
%   Example:
%
%     [I_Inn,J_Inn,A_Inn] = ...
%       assemMat_Inn_PDG2(Mesh,@STIMA_Helm_Inn_PWDG,omega);
%     [I_Imp,J_Imp,A_Imp] = ...
%       assemMat_Bnd_PDG2(Mesh,-1,@STIMA_Helm_Imp_Bnd_PWDG,omega);
%     [I_Dir,J_Dir,A_Dir] = ...
%       assemMat_Bnd_PDG2(Mesh,-2,@STIMA_Helm_Dir_Bnd_PWDG,omega);
%     A = sparse([I_Inn;I_Imp;I_Dir],...
%       [J_Inn;J_Imp;J_Dir],[A_Inn;A_Imp;A_Dir]);
%
%   See also assemMat_Bnd_PDG2, STIMA_Helm_Imp_Bnd_PWDG,
%   STIMA_Helm_Neu_Bnd_PWDG, STIMA_Helm_Inn_PWDG, LOAD_Imp_Bnd_PWDG.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize
  nDofs = Params.nDofs;             % Number of degrees of freedom
  Dir = Data.ElemData.Dir;          % Propagation directions of plane wave basis functions
  
  % Check for parametrization of edge
  isParam = isfield(Params,'Geom') && Params.Geom.isParam && nargin>=6 && ~isempty(QuadRule);
  
  % If edge is straight / not parametrized
  if(~isParam)

    % Orient normal vector as exterior normal vector
    n = -Normal;
  %   orient = check_ext_normal(Normal,Edge,Data.Vertices,Data.EdgeLoc);
  %   n = orient*Normal;

    % The fluxes will be represented as linear combinations of plane wave
    % basis functions with direction vectors Dir

    % Calculate traces (normal components)
    v = ones(nDofs,1);
    Dv = i*omega*(Dir*n.');

    % Calculate flux (flux of u is zero)
    sigma = (1/(i*omega))*Dv - Params.a*v;

    % The sought integral is the pointwise product of a factor and the L2
    % inner product matrix on the current edge.

    % Calculate integrand at first endpoint of edge
    Aloc = -i*omega*transpose(sigma(:,ones(nDofs,1))).*conj(v(:,ones(nDofs,1)));

    % Multiply by L2 inner product matrix (mass matrix)
    Aloc = Aloc.*Params.L2;
  
  % If edge is parametrized
  else
    
    % Calculate integrals directly through quadrature
    nx = numel(QuadRule.x);
    X0 = Data.Vertices(1,:);
    Gamma = Params.Geom.Gamma(QuadRule.x);
    N = Params.Geom.N(QuadRule.x);
    Aloc = Params.Geom.dGamma(ones(nDofs));
    for j=1:nDofs
      for k=1:nDofs
        ExpVal = exp(i*omega*(Gamma-X0(ones(nx,1),:))*(Dir(j,:)-Dir(k,:)).');
        Sigma = ExpVal.*(N*Dir(j,:).') - Params.a*ExpVal;
        Aloc(k,j) = -i*omega*Aloc(k,j)*sum(QuadRule.w.*Sigma);
      end
    end
    
  end
  
return