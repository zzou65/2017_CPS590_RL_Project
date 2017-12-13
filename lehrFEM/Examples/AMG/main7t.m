function [] = main7t()
% table of convergence data for discontinuous coefficients on circle

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
 
  REF = [2,6];                                   % refinements of initial mesh
  F_HANDLE = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  GD_HANDLE = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data
  M = [1,1];                                     % Smoothing steps
  
%   epsilon = 10.^(-2:-0.15:-2.45);
%   epsilon = 10.^(-1:-2:-7);
  epsilon = 10.^(-1:-1:-4);
  
  tol = 1e-6;
  maxit = 100;
  maxlev = 10;
  num = length(epsilon);
  
  AMGOpt = AMGDefaultOptions;
  AMGOpt.pre.its = M(1);
  AMGOpt.post.its = M(end);
  
  % Initialize multigrid data structure
  
  CMesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  mg_data = mg_mesh('mesh',CMesh,'ref',REF);
  mg_data = mg_smooth(mg_data,'m',M);
  mg_data = mg_error(mg_data,'energy',false,'eucl',true);
  LVL = length(mg_data);
  
  % Initialize data
  
  names = {'Geometric Multigrid','Algebraic Multigrid'};
  levels = zeros(num,2);
  dofs = zeros(maxlev,num,2);
  its = zeros(maxlev-1,num,2);
  
  for j=1:num

    % Define conductivity
    
    radius = 0.25;
    c_handle = @(x,varargin) cond(x,epsilon(j),radius);

    % Update multigrid data structure
    
    mg_data = mg_stima(mg_data,'f',F_HANDLE,'gd',GD_HANDLE,'c',c_handle);
    
    % Genereate AMG data structure
    
    AMGData = AMGSetup(mg_data{end}.A,AMGOpt);
    
    % Count degrees of freedom on each level for geometric multigrid
    
    levels(j,1) = LVL;
    for i=1:levels(j,1)
      dofs(i,j,1) = mg_data{i}.n.free;
    end
    
    % Count degrees of freedom on each level for algebraic multigrid
    
    levels(j,2) = length(AMGData);
    for i=1:levels(j,2);
      dofs(i,j,2) = size(AMGData{levels(j,2)-i+1}.A,1);
    end
    
    % calculate tolerances for geometric multigird
    % (so that geometric multigrid uses the same (weighted) norm as
    % algebraic multigird)
    
    tol_gmg = tol*4.^(LVL-2:-1:0);
    for i=1:LVL-1
      nrmb = norm(mg_data{i+1}.b);
      if(nrmb == 0)
        nrmb = 1;
      end
      tol_gmg(i) = tol_gmg(i)*nrmb;
    end
    
    % calculate solution with geometric multigrid
    
    [x,conv] = nmg_solve(mg_data,tol_gmg,maxit);
    if(conv.flag(end))
      its(1:levels(j,1)-1,j,1) = conv.iter;
    end
    
    % calculate solution with algebraic multigrid
    
    [x,flag,relres,iter] = namg_solve(mg_data{end}.b,[],AMGData,tol,maxit);
    if(flag)
      its(1:levels(j,2)-1,j,2) = iter;
    end

  end % loop over epsilons
  
  maxlev = max(levels(:));
  maxlev_gmg = max(levels(:,1));
  maxlev_amg = max(levels(:,2));
  its = its(1:maxlev-1,:,:);
  dofs = dofs(1:maxlev,:,:);
  
  % print tables
  
  fprintf('\n\n########## NUMBER OF DEGREES OF FREEDOM ##########\n\n');
  fprintf('for levels in reverse order and eps = %s\n\n',num2str(epsilon,'%10.2e'));
  
  fprintf('## %s\n\n',names{1});
  disp(dofs(maxlev_gmg:-1:1,:,1));
  
  fprintf('## %s\n\n',names{2});
  disp(dofs(maxlev_amg:-1:1,:,2));
  
  fprintf('\n\n########## NUMBER OF ITERATIONS PER LEVEL ##########\n\n');
  fprintf('for levels in reverse order and eps = %s\n\n',num2str(epsilon,'%10.2e'));
  
  fprintf('## %s\n\n',names{1});
  disp(its(maxlev_gmg-1:-1:1,:,1));

  fprintf('## %s\n\n',names{2});
  disp(its(maxlev_amg-1:-1:1,:,2));
  
  
return


% the dependence of the conductivity on epsilon is given by the function

function  c = cond(x,epsilon,r)
% construct circle with radius r

  x = x - 0.5;
  d = sqrt(sum(x.^2,2));
  ind = d < r;
  c = ones(size(x,1),1);
  c(ind) = epsilon;

return