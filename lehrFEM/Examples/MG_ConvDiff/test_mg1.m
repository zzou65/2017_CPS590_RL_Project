function [] = test_mg1()
% test geometric multigrid with linear advection

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  
  % define function handles and constants
  
  f = @(x,varargin) zeros(size(x,1),1);
  gd = @(x,varargin) double(2*x(:,1)+x(:,2)<1);
  v = @(x,varargin) ones(size(x,1),1)*[0,1];
  c = 1e-4;
  v0 = v([0 0]);
  
  tol = 1e-8;
  maxit = 50;
  
  % initialize mesh
  
  Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  Mesh = add_Edges(Mesh);
  
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Dir0 = Loc(Mesh.Coordinates(Mesh.Edges(Loc,1),2)==0 & Mesh.Coordinates(Mesh.Edges(Loc,2),2)==0);
  Dir1 = Loc(Mesh.Coordinates(Mesh.Edges(Loc,1),2)==1 & Mesh.Coordinates(Mesh.Edges(Loc,2),2)==1);
  Mesh.BdFlags(Dir0) = -1;
  Mesh.BdFlags(Dir1) = -1;
  
%   % alternatively, put Dirichlet conditions only on the inflow boundary
%   Mesh = add_Edge2Elem(Mesh);
%   [inflow,outflow,neutral] = get_inflow_outflow(Mesh,v);
%   Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
%   Mesh.BdFlags(inflow) = -1;
%   Mesh.BdFlags(outflow) = -2;
%   Mesh.BdFlags(neutral) = -2;
  
  % generate multigrid data structure
  
  mg_data = mg_mesh('mesh',Mesh,'ref',[2 6]);
  
  mg_data = mg_stima(mg_data,'stima','assem',...
    'stima_assem',@(mesh) assem_stima(mesh,c,v),...
    'f',f,'gd',gd);
  
  mg_data = mg_smooth(mg_data,'m',1,'smoother',@gs_smooth,'per','sort',...
    'per_fn',@(mesh,dofs) mesh.Coordinates(dofs,:)*v0');
  
  mg_data = mg_error(mg_data,'l2',true,'energy',false);
  
  % apply multigrid solver
  
  u = mg_data{end}.u_bd;
  [u(mg_data{end}.dofs),conv] = mg(mg_data,tol,maxit);
  
  % calculate exact solution and error
  
  u_ex = mg_data{end}.u_bd;
  u_ex(mg_data{end}.dofs) = mg_data{end}.A\mg_data{end}.b;
  err = u - u_ex;
  err_l2 = sqrt(err(mg_data{end}.dofs)'*mg_data{end}.L2*err(mg_data{end}.dofs));
  err_norm = norm(err);
  
  % print error
  
  fprintf('Error of GMG : %g (L2)  %g (Euclidean) , (%.0f iterations)\n',err_l2,err_norm,conv.iter);
  
  % plot solution

  plot_LFE(u,mg_data{end}.mesh);
  colorbar;
  title('Solution');
  
  % plot error
  
  plot_LFE(err,mg_data{end}.mesh);
  colorbar;
  set(gca,'DataAspectRatio',[1 1 max(abs(err))]);
  title('Error');
  
return 