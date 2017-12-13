% vary def. of strongly connected in coarse grid selection

%   Copyright 2006-2006 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
 
  F_HANDLE = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  GD_HANDLE = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data
  TOL = 1e-12;                                   % Stopping criterion
  MAXIT = 50;                                    % Maximum number of iterations
 
  % Generate coarse mesh

  Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  
  % Refine mesh
  
  for i = 1:5
    Mesh = refine_REG(Mesh);
  end
  
  % Construct stiffness matrix and load vector
  
  A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
  b = assemLoad_LFE(Mesh,P3O3(),F_HANDLE);
  
  % Incoporate Dirichlet boundary conditions
  
  [U,FDofs] = assemDir_LFE(Mesh,-1,GD_HANDLE);
  b = b - A*U;
  A = A(FDofs,FDofs);
  b = b(FDofs);
  
  % Define values for theta
  
  theta = [0.05,0.15,0.2,0.25,0.3,0.4,0.5];
  
  % Initialize algebraic multigrid solver
  
  n = size(theta,2);
  iter = nan(1,n);
  resvec = cell(1,n);
  time = nan(1,n);
  maxiter = 0;
  dofs = cell(1,n);
  levels = nan(1,n);
  
  opt = AMGDefaultOptions;  
  U0 = zeros(size(b));
  
  maxlvl = 0;
  
  % Run algebraic multigrid solver
  
  for i=1:n
    
    % Run amg
    
    opt.CF.theta = theta(i);
    [U(FDofs),flag,relres,iter(i),resvec{i},times,data] = amg_solve(U0,b,A,opt,TOL,MAXIT);
    time(i) = times.V_cycles.elapsed;
    maxiter = max(maxiter,iter(i));
    
    % count degrees of freedom on all levels
    
    levels(i) = length(data);
    maxlvl = max(levels(i),maxlvl);
    dofs{i} = nan(1,levels(i));
    for j=1:levels(i)
      dofs{i}(j) = size(data{j}.A,1);
    end
    
  end

  % Plot convergence
  
  figure;
  lgd = cell(1,n);
  for i=1:n
    semilogy(1:iter(i),resvec{i});
    hold all;
    lgd{i} = sprintf('%s = %3.2f (%3.2f s)','\theta',theta(i),time(i));
  end

  xlabel('\bf iteration');
  ylabel('\bf residual');
  title('\bf Convergence of Algebraic Multigrid');
  legend(lgd{:},'Location','NorthEast');
  grid on;
  set(gca,'XTick',1:maxiter,'XLim',[1,maxiter],'YLim',[1e-12,1]);
  
  
  % Plot number of degrees of freedom
  
  figure;
  lgd = cell(1,n);
  for i=1:n
    d = maxlvl-levels(i);
    plot(maxlvl:-1:d+1,dofs{i});
    hold all;
    lgd{i} = sprintf('%s = %3.2f','\theta',theta(i));
  end
  
  xlabel('\bf level');
  ylabel('\bf degrees of freedom');
  title('\bf Distribution of Degrees of Freedom for Various ''Strengths'' of Strong Connections');
  legend(lgd{:},'Location','NorthWest');
  grid on;
  set(gca,'XTick',1:maxlvl);
  
    
  % Clear memory
  
  clear all;