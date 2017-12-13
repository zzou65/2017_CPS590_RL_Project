function [] = cond_stima(ndir,refs,omega,flux_params,mtd_name)
%COND_STIMA condition number of stiffness matrix
%   
%   COND_STIME(NDIR,REFS,OMEGA,FLUX_PARAMS,MTD_NAME) plots the
%   condition number of the Helmholtz stiffness matrix for NDIR plane wave
%   basis functions, REFS mesh refinements and wave number OMEGA.  NDIR and
%   OMEGA should be vectors of the same length.
%
%   The condition number is plotted without preconditioning and with
%   preconditioning based on the singular value decompositions of the
%   energy-norm inner product and the block-diagonal of the stiffness
%   matrix.
%
%   FLUX_PARAMS is a cell array containing the flux parameter names and
%   values.  The default values are equivalent to
%       FLUX_PARAMS = {'a',0.5,'b',0.5,'c',[0 0],'d',0.5} .
%
%   MTD_NAME is a string containing the name of the method used, for
%   example 'PWDG'.
%
%   See also assemPrec_SVD_PDG.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% read arguments
if(nargin<1 || isempty(ndir))
  ndir = 3:20;
end
if(nargin<2 || isempty(refs))
  refs = 2;
end
if(nargin<3 || isempty(omega))
  omega = 2^(-1+refs)*ones(size(ndir));
elseif(numel(omega)==1)
  omega = omega(ones(size(ndir)));
end
if(nargin<4 || isempty(flux_params))
  flux_params = {};
end
if(nargin<5 || isempty(mtd_name))
  mtd_name = 'PWDG';
end

% initialize
condA = nan(size(ndir));
condApcB = nan(size(ndir));
condApcA = nan(size(ndir));
num = numel(condA);

% construct mesh
Mesh.Coordinates = [-1 -1; 1 -1; 1 1; -1 1];
Mesh.Elements = [1 2 3; 1 4 3]; % triangles
% Mesh.Elements = [1 2 3 4]; % quadrilaterals
Mesh = add_Edges(Mesh);         
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
Mesh.BdFlags(Loc) = -1;
for j=1:refs
  Mesh = refine_REG(Mesh);
end
Mesh = orient_Elems(Mesh);
Mesh = add_Edge2Elem(Mesh);
Mesh = add_DGData(Mesh);
% Mesh = set_Data_PWDG(Mesh,'Omega',omega,flux_params{:});

% loop over number of plane wave basis functions
for k=1:num
  
  % construct plane wave directions and set flux parameters
  Mesh = set_Data_PWDG(Mesh,'nDofs',ndir(k),'NewDir','replace',...
    'Omega',omega(k),flux_params{:});

  % assemble stiffness matrix
  [I_Inn,J_Inn,A_Inn] = assemMat_Inn_PDG2(Mesh,@STIMA_Helm_Inn_PWDG,omega(k));
  [I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2(Mesh,[],@STIMA_Helm_Imp_Bnd_PWDG,omega(k));
  A = sparse([I_Inn;I_Bnd],[J_Inn;J_Bnd],[A_Inn;A_Bnd]);

  % assemble energy norm inner product
  [I_inn,J_inn,B_inn] = assemMat_Inn_PDG2(Mesh,@SCAPRO_Helm_Inn_PWDG,omega(k));
  [I_bnd,J_bnd,B_bnd] = assemMat_Bnd_PDG2(Mesh,[],@SCAPRO_Helm_Bnd_PWDG,omega(k));
  B = accumarray([[I_inn;I_bnd],[J_inn;J_bnd]],[B_inn;B_bnd]);
  
  % assemble preconditioners
  [PLB,PRB] = assemPrec_SVD_PDG(Mesh,B);
  [PLA,PRA] = assemPrec_SVD_PDG(Mesh,A);
  
  % estimate condition number
  condA(k) = condest(A);            % without preconditioning
	condApcB(k) = condest(PLB*A*PRB); % with energy-norm preconditioner
  condApcA(k) = condest(PLA*A*PRA); % with stiffness-matrix preconditioner

end

% create label for omega
if(min(omega)==max(omega))
  omega_label = sprintf(', h\\omega=%g',omega(1)*2^(1-refs));
else
  omega_label = '';
end

% plot condition number
figure;
plot(ndir,condA,ndir,condApcB,ndir,condApcA);
set(gca,'YScale','log','XTick',ndir);
grid on;
legend('no prec.','energy-norm prec.','block-diag. prec.','Location','NorthWest');
xlabel('\bf number of plane wave basis functions');
ylabel('\bf condition number of stiffness matrix');
title(sprintf('\\bf Condition number of stiffness matrix for %s%s',mtd_name,omega_label));

return