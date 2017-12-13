function [] = table_linadv_bd()
% generate data tables for linear advection with boundary layer

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % define parameters

  refs = [2,6];
  CMesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  
  f = @(x,varargin) zeros(size(x,1),1);
  gd = @(x,varargin) double(2*x(:,1)+x(:,2)<1);
  v = @(x,varargin) ones(size(x,1),1)*[0,1];
  
  c = [1e-2,1e-3,1e-4,1e-5];
  
  m = 1;
  
  tol = 1e-6;
  maxit = 100;
  maxlev = 10;
  numc = length(c);
  
  AMGOpt = AMGDefaultOptions;
  AMGOpt.pre.its = m(1);
  AMGOpt.post.its = m(end);

  % define dirichlet boundary at y \in {0,1}
  
  CMesh = add_Edges(CMesh);
  Loc = get_BdEdges(CMesh);
  CMesh.BdFlags = zeros(size(CMesh.Edges,1),1);
  Dir0 = Loc(CMesh.Coordinates(CMesh.Edges(Loc,1),2)==0 ...
    & CMesh.Coordinates(CMesh.Edges(Loc,2),2)==0);
  Dir1 = Loc(CMesh.Coordinates(CMesh.Edges(Loc,1),2)==1 ...
    & CMesh.Coordinates(CMesh.Edges(Loc,2),2)==1);
  CMesh.BdFlags(Dir0) = -1;
  CMesh.BdFlags(Dir1) = -1;
  
  % initialize data
  
  names = {'Geometric Multigrid','Geometric Multigrid with Downwinding','Algebraic Multigrid','Algebraic Multigrid with Downwinding'};
  levels = zeros(numc,2);
  dofs = zeros(maxlev,numc,2);
  its = zeros(maxlev-1,numc,4);
  
  % loop over values of diffusion parameter c
  
  for k=1:numc
    
    % generate geometric multigrid data structure
  
    mg_data = mg_mesh('mesh',CMesh,'ref',refs);
    mg_data = mg_stima(mg_data,'f',f,'gd',gd,'stima','assem',...
      'stima_assem',@(mesh) assem_stima(mesh,c(k),v));
    mg_data = mg_smooth(mg_data,'m',m,'smoother',@gs_smooth);
    mg_data = mg_error(mg_data,'energy',false,'eucl',true);
    
    % generate algebraic multgrid data structure
    
    AMGData = AMGSetup(mg_data{end}.A,AMGOpt);
    
    % count degrees of freedom on each level for geometric multigrid
    
    levels(k,1) = length(mg_data);
    for i=1:levels(k,1)
      dofs(i,k,1) = mg_data{i}.n.free;
    end
    
    % count degrees of freedom on each level for algebraic multigrid
    
    levels(k,2) = length(AMGData);
    for i=1:levels(k,2);
      dofs(i,k,2) = size(AMGData{levels(k,2)-i+1}.A,1);
    end
    
    % calculate tolerances for geometric multigird
    % (so that geometric multigrid uses the same (weighted) norm as
    % algebraic multigird)
    
    LVL = length(mg_data);
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
      its(1:levels(k,1)-1,k,1) = conv.iter;
    end
    
    % calculate solution with algebraic multigrid
    
    [x,flag,relres,iter] = namg_solve(mg_data{end}.b,[],AMGData,tol,maxit);
    if(flag)
      its(1:levels(k,2)-1,k,3) = iter;
    end
    
    % construct downwind smoother for geomtric multigrid
    
    v0 = v([0 0]);
    mg_data = mg_smooth(mg_data,'per','sort',...
      'per_fn',@(mesh,dofs) mesh.Coordinates(dofs,:)*v0');
    
    % calculate solution with geometric multigrid and downwinding
    
    [x,conv] = nmg_solve(mg_data,tol_gmg,maxit);
    if(conv.flag(end))
      its(1:levels(k,1)-1,k,2) = conv.iter;
    end
    
    % construct downwind smoother for algebraic multigrid
    
    AMGOptDW = AMGOpt;
    [dummy,per] = sort(mg_data{end}.mesh.Coordinates(mg_data{end}.dofs,:)*v0');
    AMGOptDW.pre.GSperm = per;
    
    % calculate solution with algebraic multigrid and downwinding
    
    [x,flag,relres,iter] = namg_solve(mg_data{end}.b,mg_data{end}.A,AMGOptDW,tol,maxit);
    if(flag)
      its(1:levels(k,2)-1,k,4) = iter;
    end
    
    
  end % for loop over values of diffusion parameter c, index k
  
  maxlev = max(levels(:));
  maxlev_gmg = max(levels(:,1));
  maxlev_amg = max(levels(:,2));
  its = its(1:maxlev-1,:,:);
  dofs = dofs(1:maxlev,:,:);
  
  % print tables
  
  fprintf('\n\n########## NUMBER OF DEGREES OF FREEDOM ##########\n\n');
  fprintf('for levels in reverse order and c = %s\n\n',num2str(c,'%8.0e'));
  
  fprintf('## %s\n\n',names{1});
  disp(dofs(maxlev_gmg:-1:1,:,1));
  
  fprintf('## %s\n\n',names{3});
  disp(dofs(maxlev_amg:-1:1,:,2));
  
  fprintf('\n\n########## NUMBER OF ITERATIONS PER LEVEL ##########\n\n');
  fprintf('for levels in reverse order and c = %s\n\n',num2str(c,'%8.0e'));
  
  fprintf('## %s\n\n',names{1});
  disp(its(maxlev_gmg-1:-1:1,:,1));
  
  fprintf('## %s\n\n',names{2});
  disp(its(maxlev_gmg-1:-1:1,:,2));
  
  fprintf('## %s\n\n',names{3});
  disp(its(maxlev_amg-1:-1:1,:,3));
  
  fprintf('## %s\n\n',names{4});
  disp(its(maxlev_amg-1:-1:1,:,4));
  
return  