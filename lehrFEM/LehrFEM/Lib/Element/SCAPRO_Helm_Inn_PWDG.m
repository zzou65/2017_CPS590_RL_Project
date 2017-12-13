function Aloc = SCAPRO_Helm_Inn_PWDG(Edge,Normal,EdgeData,LData,RData,omega,c,varargin)
%SCAPRO_HELM_INN_PWDG Interior edge terms of PWDG inner product matrix
%
%   ALOC = SCAPRO_HELM_INN_PWDG(EDGE,NORMAL,EDGEDATA,LDATA,RDATA,OMEGA,C)
%   computes inner edge contributions to the energy norm inner product for 
%   a discontinuous plane wave discretization of the Helmholtz equation.
%
%   In the calculation of the matrix, the product of two plane waves
%   is written as a plane wave with a wave vector given by the difference
%   of the wave vectors of the plane wave basis functions (since the
%   complex conjugate of one of these is taken in the scalar product).
%   For non-diagonal entries of the matrix, this resulting plane wave
%   is then written as a constant times its Laplacian and the integral is
%   transformed into an integral of the gradient over the edge of the
%   element.  For diagonal elements, the integrand is constant and the
%   volume integral, ie. the area of the element, is calculated.
%
%   EDGE is 2-by-2 matrix whose rows contain the start and end nodes of the
%   current edge.
%
%   NORMAL is 1-by-2 marix which contains the unit normal with respect to
%   the current edge EDGE.
%
%   EDGEDATA is a structure that contains at least the fields:
%    SAMEDOFS Logical specifying whether or not the degrees of freedom on
%             the two adjacent elements are identical.
%
%   The structs LDATA and RDATA conatin the left and right hand side
%   element data:
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
%   C is a scalar factor for the L2 part of the inner product.  A value of
%   1 (default) corresponds to the energy norm, ie. the L2 part is
%   multiplied just by OMEGA squared; a value of 0 corresponds to just the
%   H1 seminorm.
%
%   Example: assembly of the energy-norm inner product matrix.
%
%     [I_inn,J_inn,B_inn] = assemMat_Inn_PDG2(Mesh,...
%       @SCAPRO_Helm_Inn_PWDG,omega,1);
%     [I_bnd,J_bnd,B_bnd] = assemMat_Bnd_PDG(Mesh,[],...
%       @SCAPRO_Helm_Bnd_PWDG,omega,1);
%     B = sparse([I_inn;I_bnd],[J_inn;J_bnd],[B_inn;B_bnd]);
%
%   See also MASS_Inn_PWDG, SCAPRO_Helm_Bnd_PWDG, assemMat_Inn_PDG2.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Set weight if it is not given
  if(nargin < 7 || isempty(c))
    c = 1;
  end
  
  % Initialize
  nDofsL = LData.ElemData.nDofs;
  DirL = LData.ElemData.Dir;
  
  nDofsR = RData.ElemData.nDofs;
  DirR = RData.ElemData.Dir;
  
  % Reduce problem to mass matrix and call corresponding code
  DirL1 = DirL(:,ones(nDofsL,1))'.*DirL(:,ones(nDofsL,1));
  DirL2 = DirL(:,2*ones(nDofsL,1))'.*DirL(:,2*ones(nDofsL,1));
  DirLMat = DirL1 + DirL2;
  
  DirR1 = DirR(:,ones(nDofsR,1))'.*DirR(:,ones(nDofsR,1));
  DirR2 = DirR(:,2*ones(nDofsR,1))'.*DirR(:,2*ones(nDofsR,1));
  DirRMat = DirR1 + DirR2;
  
  DirMat = blkdiag(DirLMat,DirRMat);

  Aloc = omega^2*(c+DirMat).*MASS_Inn_PWDG(Edge,Normal,EdgeData,LData,RData,omega,varargin{:});
  
return