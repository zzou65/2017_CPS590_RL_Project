% Test script for local multigrid
%
%   Should print:
%   levels : 8    error : 4.973038e-002    dofs : 552    iterations : 19    ssc: 6398    time : 2.26s (or whatever)
%     trunc. error : 4.589207e-003    iterations : 14    ssc : 20244    time : 0.26s (for example)

%   Copyright 2006-2006 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% initialize mesh

Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
Mesh = add_Edges(Mesh);
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = -1;
Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
for i=1:2
  Mesh = refine_REG(Mesh);
end

% run multigrid solver

[U,Mesh,flag,iter,ndof,err,timer,dofslvl,A,L,FDofs] = locmg_solve([],Mesh,@f_LShap,@g_D_LShap,'c',0.5,1,1,[1 1],@gs_smooth,0.05,50);

% print output arguments

LVL = max(Mesh.VertLevel);
fprintf('levels : %d    error : %d    dofs : %d    iterations : %d    ssc: %d    time : %4.2fs \n',...
  LVL,err(iter),ndof(iter),iter,2*sum(dofslvl(:)),timer.total);

% run multigrid solver without mesh refinements

U_bd = U;
U_bd(FDofs) = 0;
U_ex = U_bd;
U_ex(FDofs) = A(FDofs,FDofs)\(L(FDofs)-A(FDofs,:)*U_bd);
[U1,Mesh,flag1,iter1,ndof1,err1,timer1] = locmg_solve([],Mesh,L,A,U_bd,FDofs,U_ex,1,[1 1],@gs_smooth,0.005,50);

fprintf('  trunc. error : %d    iterations : %d    ssc : %d    time : %4.2fs \n',...
  err1(iter1),iter1,2*iter1*sum(dofslvl(:,end)),timer1.total);

% % plot results
% 
% plot_Mesh(Mesh,'as');
% 
% plot_LFE(U,Mesh);
% % hold on;
% % plot_Mesh(Mesh,'f');
% title('\bf Solution');
% xlabel('\bf x');
% ylabel('\bf y');
% zlabel('\bf U');
% 
% figure;
% semilogy(ndof);
% xlabel('\bf iteration');
% ylabel('\bf number of degrees of freedom N');
% title('\bf Degrees of Freedom');
% set(gca,'XLim',[1,iter]);
% grid on;
% 
% figure;
% semilogy(err);
% xlabel('\bf iteration');
% ylabel('\bf estimated error (reconstruction-based)');
% title('\bf Error');
% set(gca,'XLim',[1,iter]);
% grid on;
% 
% figure;
% its = [1,7,13,16,19];
% plot((0:LVL)',dofslvl(:,its));
% xlabel('\bf level');
% ylabel('\bf number of degrees of freedom');
% title('\bf Distribution of Degrees of Freedom');
% set(gca,'XLim',[0,LVL]);
% grid on;
% legend('iteration 1','iteration 7','iteration 13','iteration 16','iteration 19','Location','NorthEast');

% clean up

clear all;