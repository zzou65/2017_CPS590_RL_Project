function Aloc = STIMA_Helm_Inn_PWDG(Edge,Normal,Params,LData,RData,omega,varargin)
%STIMA_HELM_INN_PWDG Element stiffness matrix for PWDG Helmholtz
%
%   ALOC = STIMA_HELM_INN_PWDG(EDGE,NORMAL,PARAMS,LDATA,RDATA,OMEGA)
%   computes the entries of the element stiffness matrix for a 
%   discontinuous plane wave discretization of the Helmholtz equation on
%   interior edges.
%  
%   EDGE is 2-by-2 matrix whose rows contain the start and end nodes of the
%   current edge.
%
%   NORMAL is 1-by-2 marix which contains the unit normal with respect to
%   the current edge EDGE.
%
%   PARAMS is a structure that contains at least the fields:
%    A        Scalar coefficient for the (jump-jump) penalty term in the
%             numerical flux for the function value.
%    B        Scalar coefficient for a term containing jumps of the normal
%             derivative in the numerical flux for the gradient.
%    C        1-by-2 flux coefficient.
%    NDOFS    Total number of degrees of freedom on elements adjacent to
%             current edge.
%    L2       The L2 inner product matrix on the current edge.
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
%   Note that the interior boundary terms used here incorporate the volume
%   terms, so no separate volume terms are necessary.
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
%   See also assemMat_Inn_PDG2, STIMA_Helm_Bnd_PWDG, LOAD_Imp_Bnd_PWDG.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize
  nDofsL = LData.ElemData.nDofs;    % Number of basis functions on left element
  nDofsR = RData.ElemData.nDofs;    % Number of basis functions on right element
  nDofs = Params.nDofs;             % Total number of basis functions on left and right elements
  DirL = LData.ElemData.Dir;        % Propagation directions of basis functions on left element
  DirR = RData.ElemData.Dir;        % Propagation directions of basis functions on right element
  
  % Orient normal vector from right element to left element
  if(LData.Element>RData.Element)
    Normal = -Normal;
  end
%   orient = check_ext_normal(Normal,Edge,RData.Vertices,RData.EdgeLoc);
%   Normal = orient*Normal;

  % The fluxes will be represented as linear combinations of plane wave
  % basis functions with direction vectors Dir
  
  % Calculate traces from the left
  Lv = ones(nDofsL,1);
  LDv = i*omega*DirL;
  
  % Calculate traces from the right
  Rv = ones(nDofsR,1);
  RDv = i*omega*DirR;
  
  % Compute averages and jumps (take scalar products of vectors with normal vector)
  v = 0.5*[Lv;Rv];
  Dv = 0.5*[LDv;RDv]*Normal.';
  vN = [-Lv;Rv];
  DvN = [-LDv;RDv]*Normal.';
  cN = Params.c*Normal.';
  
  % Calculate fluxes
  sigma = 1/(i*omega)*Dv - Params.a*vN - cN/(i*omega)*DvN;
  u = v + cN*vN - Params.b/(i*omega)*DvN;
  
  % The sought integral is the pointwise product of a factor and the L2
  % inner product matrix on the current edge.
  
  % Calculate factor for integral
  Aloc = transpose(u(:,ones(nDofs,1))).*conj(DvN(:,ones(nDofs,1))) ...
    - i*omega*transpose(sigma(:,ones(nDofs,1))).*conj(vN(:,ones(nDofs,1)));

  % Multiply by L2 inner product matrix (mass matrix)
  Aloc = Aloc.*Params.L2;
  
return