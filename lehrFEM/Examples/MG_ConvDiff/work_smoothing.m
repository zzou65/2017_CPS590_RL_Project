function [] = work_smoothing()
% cost of amg v-cycle compared to gmg v-cycle

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % define parameters

  refs = [3 7];
  CMesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  
  f = @(x,varargin) zeros(size(x,1),1);
  gd = @(x,varargin) x(:,1)<0.5;
  v = @(x) ones(size(x,1),1)*[0.1,1];
  c = 1e-3;
  
  % generate multigrid data structure
  
  mg_data = mg_mesh('mesh',CMesh,'ref',refs);
  mg_data = mg_stima(mg_data,'f',f,'gd',gd,'stima','assem',...
    'stima_assem',@(mesh) assem_stima(mesh,c,v));
  
  % generate algebraic multigrid options
  
  AMGOpt = AMGDefaultOptions;
  AMGOpt.mincoarse = min(AMGOpt.mincoarse,mg_data{2}.n.free-1);

  % count degrees of freedom on each level for geometric multigrid
  
  levels_gmg = length(mg_data);
  dofs_gmg = zeros(1,levels_gmg);
  for i=1:levels_gmg
    dofs_gmg(i) = mg_data{i}.n.free;
  end
  
  % count degrees of freedom on each level for algebraic multigrid
  
  levels_amg = zeros(1,diff(refs));
  dofs_amg = cell(1,diff(refs));
  for k=1:diff(refs)
    AMGData = AMGSetup(mg_data{k+1}.A,AMGOpt);
    levels_amg(k) = length(AMGData);
    dofs_amg{k} = zeros(1,levels_amg(k));
    for i=1:levels_amg(k)
      dofs_amg{k}(i) = size(AMGData{levels_amg(k)-i+1}.A,1);
    end
  end

  % determine number of smoothing steps in a v-cycle for gmg
  
  v_gmg = zeros(1,diff(refs));
  for k=1:diff(refs)
    v_gmg(k) = sum(dofs_gmg(2:k+1));
  end
  
  % determine number of smoothing steps in a v-cycle for amg
  
  v_amg = zeros(1,diff(refs));
  for k=1:diff(refs)
    v_amg(k) = sum(dofs_amg{k}(2:levels_amg(k)));
  end
  
  % print work
  
  disp(dofs_gmg);
  
  disp(v_gmg);
  
  disp(v_amg);
  
  
return