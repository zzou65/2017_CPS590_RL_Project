function Mloc = MASS_Bnd_PWDG(Edge,Normal,EdgeData,Data,omega,QuadRule,varargin)
%MASS_BND_PWDG Boundary edge terms of element mass matrix for PWDG
%
%   MLOC = MASS_BND_PWDG(EDGE,NORMAL,EDGEDATA,DATA,OMEGA) computes the
%   contribution of boudary edges to the discontiuous plane wave mass
%   matrix.
%
%   In the calculation of the mass matrix, the product of two plane waves
%   is written as a plane wave with a wave vector given by the difference
%   of the wave vectors of the plane wave basis functions (since the
%   complex conjugate of one of these is taken in the scalar product).
%   For non-diagonal entries of the mass matrix, this resulting plane wave
%   is then written as a constant times its Laplacian and the integral is
%   transformed into an integral of the gradient over the edge of the
%   element.  For diagonal elements, the integrand is constant and the
%   volume integral, ie. the area of the element, is calculated.
%
%   EDGE is 2-by-2 matrix whose rows contain the start and end node of the
%   current edge.
%
%   NORMAL is 1-by-2 marix which contains the interior unit normal with
%   respect to the current edge EDGE.
%
%   EDGEDATA is a structure that contains at least the fields:
%    L2       The L2 inner product matrix on the current edge.
%    LENGTH   Length of the current edge.
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
%   DIR is an N-by-2 matrix containing the direction vectors of the plane
%   wave basis functions.
%
%   Example: assembly of the mass matrix.
%
%     [I_inn,J_inn,M_inn] = assemMat_Inn_PDG2(Mesh,...
%       @MASS_Inn_PWDG,omega);
%     [I_bnd,J_bnd,M_bnd] = assemMat_Bnd_PDG2(Mesh,[],...
%       @MASS_Bnd_PWDG,omega);
%     M = sparse([I_inn;I_bnd],[J_inn;J_bnd],[M_inn;M_bnd]);
%
%   See also MASS_Inn_PWDG, assemMat_Bnd_PDG2.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize
  nDofs = Data.ElemData.nDofs;      % Number of degrees of freedom on adjacent element
  Dir = Data.ElemData.Dir;          % Propagation directions of basis functions
  Mloc = zeros(nDofs);
  
  % Check for parametrization of edge
  isParam = isfield(EdgeData,'Geom') && EdgeData.Geom.isParam && nargin>=6 && ~isempty(QuadRule);
  
  % If edge is straight / not parametrized
  if(~isParam)
  
    % Construct indices for off-diagonal and strictly upper triangular
    % entries
    I = true(nDofs);
    upper = triu(I,1);
    offdiag = upper | upper.';

    % Calculate vector tangent to the edge and its norm
    BmA = Edge(2,:) - Edge(1,:);
    nrmBmA = EdgeData.Length;

    % Project the propagation directions of plane wave basis functions onto
    % the edge and normal vectors
    DirN = Dir*Normal';
    dDirN = DirN(:,ones(nDofs,1))' - DirN(:,ones(nDofs,1));

    DirBmA = Dir*BmA';
    dDirBmA = DirBmA(:,ones(nDofs,1))' - DirBmA(:,ones(nDofs,1));

    % Compute norm squared of propagation directions
    nrmdDir2 = dDirN.^2 + dDirBmA.^2/nrmBmA^2;

    % Construct constant factor for mass matrix
    Mloc(offdiag) = (i/omega)*(dDirN(offdiag)./nrmdDir2(offdiag));

    % Multiply by L2 inner product on edge
    Mloc = Mloc.*EdgeData.L2;
    
    % Compute diagonal entries
    Mloc(~offdiag) = -nrmBmA*(Edge(1,1) + 0.5*BmA(1))*Normal(1);
  
  % If edge is parametrized
  else
    
    % Calculate integrals directly through quadrature
    nx = numel(QuadRule.x);
    X0 = Data.Vertices(1,:);
    Gamma = EdgeData.Geom.Gamma(QuadRule.x);
    N = EdgeData.Geom.N(QuadRule.x);
    D = EdgeData.Geom.dGamma*sum(QuadRule.w.*Gamma(:,1).*N(:,1));
    for j=1:nDofs
      Mloc(j,j) = D;
      for k=j+1:nDofs
        d = Dir(j,:)-Dir(k,:);
        ExpVal = exp(i*omega*(Gamma-X0(ones(nx,1),:))*d.');
        Mloc(k,j) = EdgeData.Geom.dGamma/(i*omega*sum(d.^2))...
          *sum(QuadRule.w.*(N*d.').*ExpVal);
        Mloc(j,k) = conj(Mloc(k,j));
      end
    end
    
  end
  
return