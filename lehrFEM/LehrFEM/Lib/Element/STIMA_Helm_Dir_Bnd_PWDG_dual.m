function Aloc = STIMA_Helm_Dir_Bnd_PWDG_dual(Edge,Normal,Params,varargin)
%STIMA_HELM_DIR_BND_PWDG_DUAL Element stiffness matrix for mixed PWDG
% 
%   ALOC = STIMA_HELM_DIR_BND_PWDG_DUAL(EDGE,NORMAL,PARAMS)
%   computes the Dirichlet boundary edge contributions to the element
%   stiffness matrix for a dual (mixed) discontinuous plane wave
%   discretization of the  Helmholz equation.
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
  nDofs = Params.nDofs;             % Number of degrees of freedom
  
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
  v = [ones(nDofs,1);zeros(2*nDofs,1)];
  tau = [zeros(nDofs,1);ones(2*nDofs,1)].*N;
  
  % Calculate flux (flux of u is zero)
  sigma = tau - Params.a*v;
  
  % The sought integral is the pointwise product of a factor and the L2
  % inner product matrix on the current edge.
  
  % Calculate factor for integral
  Aloc = -transpose(sigma(:,ones(3*nDofs,1))).*conj(v(:,ones(3*nDofs,1)));
  
  % Construct integral over edge from local L2 inner product
  L2 = repmat(Params.L2,3,3);
  
  % Multiply to get edge contribution to stiffness matrix
  Aloc = Aloc.*L2;

return