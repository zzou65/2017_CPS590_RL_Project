function [] = maina()
% convergence of amg on refined meshes
%
%   PROBLEM : the vector iteration does NOT converge; eigenvalues of error
%   propagation matrix are probably incorrect !!!!!!!!!!!!!!!!!!!!!!!!!!!!

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize parameters
  
  refs = 0:5:100;          % number of local refinements
  drefs = diff([0,refs]);
  maxit = 30;               % number of iterations in power method

  % Initialize mesh
  
  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);
  for i=1:4
    Mesh = refine_REG(Mesh);
  end
  Mesh = init_LEB(Mesh);
  Mesh = rmfield(Mesh,'BdFlags');
  
  % Initialize solution
  
  c_amg = nan(size(refs));
  levels = nan(size(refs));
  
  for j=1:numel(refs)
  
    % Do local refinements towards 0

    for i=1:drefs(j)
      Mesh = add_Patches(Mesh);
      eta = pointref([],Mesh);
      marked = find(eta);
      Mesh = refine_LEB(Mesh,marked);
      Mesh = add_Edges(Mesh);
    end


    % Assemble stiffness matrix

    A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);


    % Generate AMG data structure

    AMGOptions = AMGDefaultOptions;
    AMGOptions.mincoarse = 10;
    
    AMGData = AMGSetup(A,AMGOptions);

    levels(j) = length(AMGData);
    

    % Calculate convergence rate for algebraic multigrid using the power
    % method
    %
    % DOES NOT CONVERGE !!!

    Eerr = ones(size(A,1),1);
    for i=1:maxit
      err = Eerr;
      c = AMGVcycle(AMGData,A*err);
      Eerr = err - c;
      c_amg(j) = (err'*Eerr)/(err'*err);
    end
  
  end % for loop with index j
  
  
  % Plot maximal eigenvalue of error propagation matrix against number of
  % local refinements
  
  figure;
  plot(refs,c_amg,'-+');
  xlabel('\bf number of refinements towards 0');
  ylabel('\bf maximal eigenvalue of error propagation matrix');
  title('\bf Convergence of AMG for Local Refinements');
  grid on;
  
  % Plot number of levels in AMG against the number of local refinements
  
  figure;
  plot(refs,levels,'-+');
  xlabel('\bf number of refinements towards 0');
  ylabel('\bf number of levels');
  title('\bf Number of Levels in AMG for Local Refinements');
  grid on;
  set(gca,'YLim',[1,max(levels+1)]);

return