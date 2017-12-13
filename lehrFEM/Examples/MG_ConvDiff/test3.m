function [] = test3()
% test algebraic multigrid with linear advection

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % define parameters

  refs = [2,6];
  Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  
  handle.f = @(x,varargin) zeros(size(x,1),1);
%   handle.gd = @(x,varargin) x(:,1)<0.5;
%   handle.v = @(x) ones(size(x,1),1)*[0.1,1];
  handle.gd = @(x,varargin) double(2*x(:,1)+x(:,2)<1);
  handle.v = @(x,varargin) ones(size(x,1),1)*[0,1];
  handle.c = 1e-4;
  
  m = 1;
  tol = 1e-6;
  maxit = 50;

  % initialize and refine mesh data structure
  
  Mesh = add_Edges(Mesh);
  Mesh = add_Edge2Elem(Mesh);
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
%   Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  
  % define dirichlet boundary at y \in {0,1}
  
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Dir0 = Loc(Mesh.Coordinates(Mesh.Edges(Loc,1),2)==0 & Mesh.Coordinates(Mesh.Edges(Loc,2),2)==0);
  Dir1 = Loc(Mesh.Coordinates(Mesh.Edges(Loc,1),2)==1 & Mesh.Coordinates(Mesh.Edges(Loc,2),2)==1);
  Mesh.BdFlags(Dir0) = -1;
  Mesh.BdFlags(Dir1) = -1;
  
  for k=1:refs(1)
    Mesh = refine_REG(Mesh);
  end
  
%   Mesh = rmfield(Mesh,'BdFlags');
%   Mesh = add_Edge2Elem(Mesh);
%   [inflow,outflow,neutral] = get_inflow_outflow(Mesh,handle.v);
%   Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
%   Mesh.BdFlags(inflow) = -1;
%   Mesh.BdFlags(outflow) = -2;
%   Mesh.BdFlags(neutral) = -2;
  
  for k=refs(1)+1:refs(2)
    Mesh = refine_REG(Mesh);
  end
  
  % compute stiffness matrix on finest level
  
  if(isa(handle.c,'numeric'))
    A = handle.c*assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
  else
    A = assemMat_LFE(Mesh,@STIMA_Heat_LFE,P1O2(),handle.c);
  end
  M = assemMat_MassZeroD(Mesh);
  D_fd = assemMat_LFE(Mesh,@LIEUP_FD,handle.v);

  A_fd = A+M*D_fd;
  
  % compute load vector and incorporate Dirichlet boundary
  
  L = assemLoad_LFE(Mesh,P7O6(),handle.f);
  [U_bd,FDofs] = assemDir_LFE(Mesh,-1,handle.gd);
  L = L-A_fd*U_bd;
  
  % define matrix and right hand side on free degrees of freedom
  
  iA = A_fd(FDofs,FDofs);
  iL = L(FDofs);
  
  % get exact solution
  
  U_ex = U_bd;
  U_ex(FDofs) = iA\iL;
  
  % define AMG options
  
  AMGOpt = AMGDefaultOptions;
  AMGOpt.pre.its = m(1);
  AMGOpt.post.its = m(end);
%   AMGOpt.levsmax = 2;

  v = handle.v([0,0]);
  [dummy,per] = sort(v(1)*Mesh.Coordinates(FDofs,1)+v(2)*Mesh.Coordinates(FDofs,2));
  AMGOpt.pre.GSperm=per;
  
  % generate AMG data
  
  AMGData = AMGSetup(iA,AMGOpt);
  
  % solve sample problem
  
  U = U_bd;
  [U(FDofs),flag,relres,iter] = amg_solve(U(FDofs),iL,[],AMGData,tol,maxit);
  
  err = norm(U-U_ex);
  fprintf('Error of AMG : %g , (%.0f iterations)\n',err,iter);
  
  % plot solution

  plot_LFE(U,Mesh);
  colorbar;
  
  % plot error
  
  Err = U-U_ex;
  plot_LFE(Err,Mesh);
  colorbar;
  set(gca,'DataAspectRatio',[1 1 max(abs(Err))]);
  
  % plot coarse grid points
  
  plot_AMG_coarse(Mesh,AMGData,FDofs);

return
  