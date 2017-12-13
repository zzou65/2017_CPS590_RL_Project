function Aloc = STIMA_Helm_Vol_PWDG_dual(Vertices,ElemData,omega,M_)
%STIMA_HELM_VOL_PWDG_DUAL Element stiffness matrix for mixed PWDG
%
%   ALOC = STIMA_HELM_VOL_PWDG_DUAL(VERTICES,ELEMDATA,OMEGA,M)
%   computes the volume contributions the element stiffness matrix for a
%   dual (mixed) discontinuous plane wave discretization of the Helmholz
%   equation.
%
%   VERTICES is 3-by-2 or 4-by-2 matrix specifying the vertices of the
%   current element in a row-wise fashion.
%
%   ELEMDATA is a structure that contains at least the fields:
%     NDOFS   The number of degrees of freedom on the current element.
%     DIR     A NDOFS-by-2 matrix containing the propagation directions of
%             the plane wave basis functions in its rows.
%     IND     Vector of length NDOFS containing the global indices of the
%             local degrees of freedom.
%
%   OMEGA is the wave number of the Helholtz equation.
%
%   M is the global mass matrix.
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
  
  % The functions will be represented as linear combinations of plane wave
  % basis functions with direction vectors Dir.  Note that, since the
  % values are three-dimensional vectors (function & gradient), this
  % representation will have dimension 3*nDofs.

  % Initialize
  nDofs = ElemData.nDofs;           % Number of degrees of freedom
  Dir = ElemData.Dir;               % Propagation directions of plane wave basis functions
        
  % Construct divergence and gradient
  divTau1 = -i*omega*Dir(:,1*ones(1,nDofs));
  divTau2 = -i*omega*Dir(:,2*ones(1,nDofs));
  gradV1 = -i*omega*Dir(:,1*ones(1,nDofs));
  gradV2 = -i*omega*Dir(:,2*ones(1,nDofs));
  c = i*omega(ones(nDofs));
  o = zeros(nDofs);
  
  % Construct integrand
  Aloc = [c gradV1 gradV2
          divTau1 c o
          divTau2 o c ];
        
  % Construct integral over element from mass matrix
  M0 = M_(ElemData.Ind,ElemData.Ind);
  M = repmat(M0,3,3);
  
  % Multiply to get volume contribution to stiffness matrix
  Aloc = Aloc.*M;
  
return