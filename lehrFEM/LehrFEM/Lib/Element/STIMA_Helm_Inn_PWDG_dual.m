function Aloc = STIMA_Helm_Inn_PWDG_dual(Edge,Normal,Params,LData,RData,varargin)
%STIMA_HELM_INN_PWDG_DUAL Element stiffness matrix for mixed PWDG
%
%   ALOC = STIMA_HELM_INN_PWDG_DUAL(EDGE,NORMAL,PARAMS,LDATA,RDATA)
%   computes the interior edge contributions the element stiffness matrix
%   for a dual (mixed) discontinuous plane wave discretization of the
%   Helmholz equation.
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
%   Example:
%
%     [I_Vol,J_Vol,A_Vol] = assemMat_Vol_PDG2_vec(...
%       Mesh,3,@STIMA_Helm_Vol_PWDG_dual,omega,M);
%     [I_Inn,J_Inn,A_Inn] = assemMat_Inn_PDG2_vec(...
%       Mesh,3,@STIMA_Helm_Inn_PWDG_dual);
%     [I_Imp,J_Imp,A_Imp] = assemMat_Bnd_PDG2_vec(...
%       Mesh,-1,3,@STIMA_Helm_Imp_Bnd_PWDG_dual);
%     [I_Dir,J_Dir,A_Dir] = assemMat_Bnd_PDG2_vec(...
%       Mesh,-2,3,@STIMA_Helm_Dir_Bnd_PWDG_dual);
%     A = sparse([I_Vol;I_Inn;I_Imp;I_Dir],...
%       [J_Vol;J_Inn;J_Imp;J_Dir],[A_Vol;A_Inn;A_Imp;A_Dir]);

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize
  nDofsL = LData.ElemData.nDofs;    % Number of basis functions on left element
  nDofsR = RData.ElemData.nDofs;    % Number of basis functions on right element
  nDofs = Params.nDofs;             % Total number of basis functions on left and right elements
  
  % Orient normal vector from right element to left element
  if(LData.Element>RData.Element)
    Normal = -Normal;
  end
%   orient = check_ext_normal(Normal,Edge,RData.Vertices,RData.EdgeLoc);
%   Normal = orient*Normal;
  
  % The fluxes will be represented as linear combinations of plane wave
  % basis functions with direction vectors [DirL;DirR].  Note that, since
  % the values are three-dimensional vectors (function & gradient), this
  % representation will have dimension 3*nDofs.

  % Constuct normal vector and parameter c in vector form
  N = [zeros(nDofs,1);Normal(ones(nDofs,1),1);Normal(ones(nDofs,1),2)];
  cN = Params.c*Normal.';
  
  % Calculate traces from left
  Lv = [ones(nDofsL,1);zeros(2*nDofsL+3*nDofsR,1)];
  Ltau = [zeros(nDofsL+nDofsR,1);ones(nDofsL,1);zeros(nDofsR,1);ones(nDofsL,1);zeros(nDofsR,1)];
  
  % Calculate traces from right
  Rv = [zeros(nDofsL,1);ones(nDofsR,1);zeros(2*nDofsL+2*nDofsR,1)];
  Rtau = [zeros(2*nDofsL+nDofsR,1);ones(nDofsR,1);zeros(nDofsL,1);ones(nDofsR,1)];
  
  % Compute averages and jumps (take scalar products of vectors with normal vector)
  v = 0.5*(Lv + Rv);
  tau = 0.5*(Ltau + Rtau).*N;
  v_N = (-Lv + Rv);
  tau_N = (-Ltau + Rtau).*N;

  % Calculate fluxes
  u = v + cN*v_N - Params.b*tau_N;
  sigma = tau - Params.a*v_N - cN*tau_N;
  
  % The sought integral is the pointwise product of a factor and the L2
  % inner product matrix on the current edge.
  
  % Calculate factor for integral
  Aloc1 = transpose(u(:,ones(3*nDofs,1))).*conj(tau_N(:,ones(3*nDofs,1)));
  Aloc2 = transpose(sigma(:,ones(3*nDofs,1))).*conj(v_N(:,ones(3*nDofs,1)));
  Aloc = -Aloc1-Aloc2;
  
  % Construct integral over edge from local L2 inner product
  L2 = repmat(Params.L2,3,3);
  
  % Multiply to get edge contribution to stiffness matrix
  Aloc = Aloc.*L2;

return