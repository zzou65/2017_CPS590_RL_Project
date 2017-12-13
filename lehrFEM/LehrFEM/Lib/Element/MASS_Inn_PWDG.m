function Mloc = MASS_Inn_PWDG(Edge,Normal,EdgeData,LData,RData,omega,varargin)
%MASS_INN_PWDG Interior edge terms of element mass matrix for PWDG
%
%   MLOC = MASS_INN_PWDG(EDGE,NORMAL,EDGEDATA,LDATA,RDATA,OMEGA) computes
%   the contribution of interior edges to the discontiuous plane wave mass
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
%   EDGE is 2-by-2 matrix whose rows contain the start and end nodes of the
%   current edge.
%
%   NORMAL is 1-by-2 marix which contains the unit normal with respect to
%   the current edge EDGE.
%
%   EDGEDATA is a structure that contains at least the fields:
%    SAMEDOFS Logical specifying whether or not the degrees of freedom on
%             the two adjacent elements are identical.
%    L2       The L2 inner product matrix on the current edge.
%    LENGTH   Length of the current edge.
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
%   Example: assembly of the mass matrix.
%
%     [I_inn,J_inn,M_inn] = assemMat_Inn_PDG2(Mesh,...
%       @MASS_Inn_PWDG,omega);
%     [I_bnd,J_bnd,M_bnd] = assemMat_Bnd_PDG2(Mesh,[],...
%       @MASS_Bnd_PWDG,omega);
%     M = sparse([I_inn;I_bnd],[J_inn;J_bnd],[M_inn;M_bnd]);
%
%   See also MASS_Bnd_PWDG, assemMat_Inn_PDG2.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize
  nDofsL = LData.ElemData.nDofs;    % Number of degrees of freedom on left element
  DirL = LData.ElemData.Dir;        % Propagation directions of basis functions on left element

  nDofsR = RData.ElemData.nDofs;    % Number of degrees of freedom on right element
  DirR = RData.ElemData.Dir;        % Propagation directions of basis functions on right element
  doR = ~EdgeData.sameDofs;         % If the dofs are identical, some computations become redundant
  
  % Preallocate memory
  MLint = zeros(nDofsL); 
  MRint = zeros(nDofsR);
  
  % Construct indices for off-diagonal and strictly upper triangular
  % entries for left and right elements
  I = true(nDofsL);
  upperL = triu(I,1);
  offdiagL = upperL | upperL.';
  
  I = true(nDofsR);
  upperR = triu(I,1);
  offdiagR = upperR | upperR.';  
  
  % Orient normal vector from right element to left element
  if(LData.Element>RData.Element)
    Normal = -Normal;
  end
  
  % Calculate vector tangent to the edge and its norm
  BmA = Edge(2,:) - Edge(1,:);
  nrmBmA = EdgeData.Length;
  
  % Project the propagation directions of plane wave basis functions onto
  % the edge and normal vectors for the left element
  DirLN = DirL*Normal';
  dDirLN = DirLN(:,ones(nDofsL,1))' - DirLN(:,ones(nDofsL,1));
  
  DirLBmA = DirL*BmA';
  dDirLBmA = DirLBmA(:,ones(nDofsL,1))' - DirLBmA(:,ones(nDofsL,1));
    
  % Compute norm squared of propagation directions for left element
  nrmdDirL2 = dDirLN.^2 + dDirLBmA.^2/nrmBmA^2;
   
  % Construct constant factor for mass matrix on left element
  MLint(offdiagL) = (i/omega)*(dDirLN(offdiagL)./nrmdDirL2(offdiagL));
  
  % Repeat the above for the right element
  if(doR)
    DirRN = DirR*Normal';
    dDirRN = DirRN(:,ones(nDofsR,1))' - DirRN(:,ones(nDofsR,1));

    DirRBmA = DirR*BmA';
    dDirRBmA = DirRBmA(:,ones(nDofsR,1))' - DirRBmA(:,ones(nDofsR,1));
    
    nrmdDirR2 = dDirRN.^2 + dDirRBmA.^2/nrmBmA^2;
  
    MRint(offdiagR) = (-i/omega)*(dDirRN(offdiagR)./nrmdDirR2(offdiagR));
  else
    MRint = -MLint;
  end
  
  % Construct local mass matrix
  Mloc = zeros(nDofsL+nDofsR);
  indL = 1:nDofsL;
  indR = nDofsL+1:nDofsL+nDofsR;
  Mloc(indL,indL) = MLint;
  Mloc(indR,indR) = MRint;
  Mloc = Mloc.*EdgeData.L2;
  
  % Compute contribution to diagonal entries
  d = nrmBmA*(Edge(1,1) + 0.5*BmA(1))*Normal(1);
  Mloc(indL,indL) = Mloc(indL,indL) - diag(d(ones(nDofsL,1)));
  Mloc(indR,indR) = Mloc(indR,indR) + diag(d(ones(nDofsR,1)));

return