function [] = maind()
% cost of amg v-cycle compared to gmg v-cycle

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
 
  REF = [3,7];                                   % refinements of initial mesh
  F_HANDLE = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  GD_HANDLE = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data
  
  % Generate multigrid data structure
  
  CMesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  mg_data = mg_mesh('mesh',CMesh,'ref',REF);
  mg_data = mg_stima(mg_data,'f',F_HANDLE,'gd',GD_HANDLE);
  
  % Generate algebraic multigrid options
  
  AMGOpt = AMGDefaultOptions;
  AMGOpt.mincoarse = min(AMGOpt.mincoarse,mg_data{2}.n.free-1);
  
  % Count degrees of freedom on each level for geometric multigrid
  
  levels_gmg = length(mg_data);
  dofs_gmg = zeros(1,levels_gmg);
  for i=1:levels_gmg
    dofs_gmg(i) = mg_data{i}.n.free;
  end
  
  % Count degrees of freedom on each level for algebraic multigrid
  
  levels_amg = zeros(1,diff(REF));
  dofs_amg = cell(1,diff(REF));
  for k=1:diff(REF)
    AMGData = AMGSetup(mg_data{k+1}.A,AMGOpt);
    levels_amg(k) = length(AMGData);
    dofs_amg{k} = zeros(1,levels_amg(k));
    for i=1:levels_amg(k)
      dofs_amg{k}(i) = size(AMGData{levels_amg(k)-i+1}.A,1);
    end
  end
  
%   % Print number of degrees of freedom
%   
%   dofs_gmg
%   
%   for k=1:diff(REF)
%     dofs_amg{k}
%   end
%   
  % Determine number of smoothing steps in a v-cycle for gmg
  
  v_gmg = zeros(1,diff(REF));
  for k=1:diff(REF)
    v_gmg(k) = sum(dofs_gmg(2:k+1));
  end
  
  % Determine number of smoothing steps in a v-cycle for amg
  
  v_amg = zeros(1,diff(REF));
  for k=1:diff(REF)
    v_amg(k) = sum(dofs_amg{k}(2:levels_amg(k)));
  end
  
  % Print work
  
  disp(dofs_gmg);
  
  disp(v_gmg);
  
  disp(v_amg);
  
  
return