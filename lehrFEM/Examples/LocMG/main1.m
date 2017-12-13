% Convergence of local refinement using multigrid
%
%   The effect of using multigrid is that only a limited amount of work is
%   done before a refinement.

%   Copyright 2006-2006 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% initialize coarse mesh

Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;
Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
for i=1:2
  Mesh = refine_REG(Mesh);
end

% run local multigrid solver

theta = 0.5;
[U,FMesh,flag,iter,ndof,err] = locmg_solve([],Mesh,@f_LShap,@g_D_LShap,'c',theta,1,1,[1,1],@gs_smooth,0.02,50);

% run normal multigrid solver

maxit = iter;
lvl = max(FMesh.VertLevel);
itlvl = ceil(maxit/min(lvl,3));
[U_mg,FMesh_mg,flag_mg,iter_mg,ndof_mg,err_mg] = locmg_solve([],Mesh,@f_LShap,@g_D_LShap,'c',0,itlvl,1,[1,1],@gs_smooth,0.07,maxit);

% plot error vs degrees of freedom

figure;
loglog(ndof,err,ndof_mg,err_mg);

p = polyfit(log(ndof),log(err),1);
add_Slope(gca,'SouthWest',p(1));

p = polyfit(log(ndof_mg),log(err_mg),1);
add_Slope(gca,'NorthWest',p(1));

xlabel('\bf N');
ylabel('\bf estimeted error (reconstruction-based)');
title('\bf Convergence of Local Refinements');
legend('Local Multigrid','Multigrid','Location','NorthEast');
grid on;

% clean up

clear all;
