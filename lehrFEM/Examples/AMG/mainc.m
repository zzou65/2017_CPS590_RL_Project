function [] = mainc()
% fillup of stiffness matrices on coarse grids
%
%   Determines the average number of non-zero elements per row of the
%   stiffness matrices on each level of geometric multigrid and algebraic
%   multigrid.  For AMG, the (Ruge-Stueben) interpolation is done 
%   on all coarse grid points and separately using only coarse grid points
%   to which a given fine grid point is strongly connected.  In the first
%   case, there is a fill-up in the stiffness matrix; the second version of
%   AMG has a constant number of non-zero entries per row which is about
%   twice that of geometric multigrid.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

 
  % Generate coarse mesh

  Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  
  % Initialize geometric multigrid
  
  f_handle = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  gd_handle = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data

  mg_data = mg_mesh('mesh',Mesh,'ref',[2 7]);
  mg_data = mg_stima(mg_data,'f',f_handle,'gd',gd_handle);
  
  % Generate AMG data structures
  
  A = mg_data{end}.A;
  
  AMGOptions = AMGDefaultOptions;
  AMGData1 = AMGSetup(A,AMGOptions);
  
  AMGOptions.P.nofill64 = 0;
  AMGData2 = AMGSetup(A,AMGOptions);

  % Determine the average number of non-zero entries per line for geometric
  % multigrid
  
  LVL_mg = length(mg_data);
  nz_mg = nan(1,LVL_mg);
  for j=1:LVL_mg
    nz_mg(j) = nnz(mg_data{LVL_mg+1-j}.A)/mg_data{LVL_mg+1-j}.n.free;
  end
  
  % Determine the average number of non-zero entries per line for algebraic
  % multigrid with no control of fillup
  
  LVL_amg = length(AMGData2);
  nz_amg = nan(1,LVL_amg);
  for j=1:LVL_amg
    nz_amg(j) = nnz(AMGData2{j}.A)/size(AMGData2{j}.A,1);
  end
  
  % Determine the average number of non-zero entries per line for algebraic
  % multigrid using only strong connections in interpolation
  
  LVL_amg_nofill = length(AMGData1);
  nz_amg_nofill = nan(1,LVL_amg_nofill);
  for j=1:LVL_amg_nofill
    nz_amg_nofill(j) = nnz(AMGData1{j}.A)/size(AMGData1{j}.A,1);
  end
  
  % Plot average number of non-zero entries per line of stiffness matrix
  % against the level for all three methods
  
  figure;
  plot(LVL_mg:-1:1,nz_mg,'-o',...
       LVL_amg:-1:1,nz_amg,'-+',...
       LVL_amg_nofill:-1:1,nz_amg_nofill,'-x');
  set(gca,'XDir','reverse','XTick',1:max([LVL_mg,LVL_amg,LVL_amg_nofill]));
  grid on;
  title('\bf Fill-Up of Stiffness Matrix in AMG');
  xlabel('\bf level');
  ylabel('\bf average number of non-zero entries per row');
  legend('Geometric MG','AMG (all con.)','AMG (strong con.)','Location','NorthWest');
  
return