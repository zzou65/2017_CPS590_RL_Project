% Run script for several refinement strategies

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize constants
  
  F_HANDLE = @f_LShap;      % Right hand side source term
  GD_HANDLE = @g_D_LShap;   % Dirichlet boundary data
  UEX_HANDLE = @uex_LShap;  % Exact siolution
  TOL = .01;                % Stopping criterion
  THETA = .5;               % Refinement strategy
  NITER = 5;                % Plotting levels
  
  % Initialize mesh
  
  Mesh=load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
  Mesh.BdFlags(Loc) = -1;
  Mesh = init_LEB(Mesh);
      
  iter = 0;
  Eta = TOL+1;
  d=size(Mesh.Coordinates,1);
  nmax=0;
  while(d < nmax)
    
 %   plot_Mesh(Mesh,'petas');  
    iter = iter+1;
    nDofs(iter) = size(Mesh.Coordinates,1);
    d=nDofs(iter);
      
    % Compute FE solution
  
    A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
    L = assemLoad_LFE(Mesh,P7O6(),F_HANDLE);
    [U,FreeDofs] = assemDir_LFE(Mesh,-1,GD_HANDLE);
    L = L - A*U;
    if(~isempty(FreeDofs))
      U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    end
    
    % Real and estimated error
    QuadRule=Duffy(TProd(gauleg(1,0,10)));
    Err(iter) = H1SErr_LFE(Mesh,U,QuadRule,UEX_HANDLE);  
    Mesh = add_Patches(Mesh);
    Eta_K = H1SErrDistr_LFE(Mesh,U,QuadRule,UEX_HANDLE);
    
    
    Eta(iter) = sqrt(sum(Eta_K.^2));
    %    plot_EST(Eta_K,Mesh);
    % Mark elements for refinement
    
    % first marking strategy
    Eta_max = max(Eta_K);
    MarkedElem = find(Eta_K > Eta_max*THETA);
%     
    % second marking strategy
%     [Eta_K_S,id]=sort(Eta_K,'descend');
%     csum=sqrt(sum(Eta_K.^2));
%     psum=Eta_K_S(1);
%     i=1;
%     while (psum<THETA*csum)
%         i=i+1;
%         psum=psum+Eta_K_S(i);
%     end;    
%     MarkedElem=id(1:i);    
    % Refine mesh by largest edge bisection
    
    tmp_Mesh = Mesh;
    Mesh = refine_LEB(Mesh,MarkedElem);
    
    % Update mesh data structure
    
    Mesh = add_Edges(Mesh);
    Mesh = rmfield(Mesh,'BdFlags');
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
    Mesh.BdFlags(Loc) = -1;   
  
    % Print out information
    
    fprintf('Eta  :  %2.4e,  nDofs  :  %2.4e\n',Err(iter),d);
     % Print out mesh if not to big
    if (nDofs<750)
     figure;
     plot_Mesh(Mesh,'as');
     file = ['Mesh_max_' int2str(nDofs(iter)) '.eps'];
     print('-depsc',file);
    end
  end
  
    % second marking strategy
    % Initialize mesh
  
  Mesh=load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
  Mesh.BdFlags(Loc) = -1;
  Mesh = init_LEB(Mesh);
      
  
  iter1 = 0;
  Eta = TOL+1;
  d=size(Mesh.Coordinates,1);
  dmax=0;
  while(d < dmax)
    
 %   plot_Mesh(Mesh,'petas');  
    iter1 = iter1+1;
    nDofs1(iter1) = size(Mesh.Coordinates,1);  
    d=nDofs1(iter1)  ;
    % Compute FE solution
  
    A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
    L = assemLoad_LFE(Mesh,P7O6(),F_HANDLE);
    [U,FreeDofs] = assemDir_LFE(Mesh,-1,GD_HANDLE);
    L = L - A*U;
    if(~isempty(FreeDofs))
      U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    end
    
    % Real and estimated error
    QuadRule=Duffy(TProd(gauleg(1,0,10)));
    Err1(iter1) = H1SErr_LFE(Mesh,U,QuadRule,UEX_HANDLE);  
    Mesh = add_Patches(Mesh);
    Eta_K = H1SErrDistr_LFE(Mesh,U,QuadRule,UEX_HANDLE);
    
    % Mark elements for refinement
    
    % first marking strategy
%     Eta_max = max(Eta_K);
%     MarkedElem = find(Eta_K > Eta_max*THETA);
%     
    % second marking strategy
    [Eta_K_S,id]=sort(Eta_K,'descend');
    csum=sqrt(sum(Eta_K.^2));
    psum=Eta_K_S(1);
    i=1;
    while (psum<THETA*csum)
        i=i+1;
        psum=psum+Eta_K_S(i);
    end;    
    MarkedElem=id(1:i);    
    % Refine mesh by largest edge bisection
    
    tmp_Mesh = Mesh;
    Mesh = refine_LEB(Mesh,MarkedElem);
    
    % Update mesh data structure
    
    Mesh = add_Edges(Mesh);
    Mesh = rmfield(Mesh,'BdFlags');
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
    Mesh.BdFlags(Loc) = -1;   
  
    % Print out information
    
    fprintf('Eta  :  %2.4e,  nDofs  :  %2.4e\n',Err1(iter1),d);
    % Print out mesh if not to big
    if (nDofs1<750)
     figure;
     plot_Mesh(Mesh,'as');
     file = ['Mesh_sum_' int2str(nDofs1(iter1)) '.eps'];
     print('-depsc',file);
    end
  end
  
 % regularrefinement
 % Initialize mesh
  
  Mesh=load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
  Mesh.BdFlags(Loc) = -1;
  Mesh = init_LEB(Mesh);
      
  
  iter2 = 0;
  Eta = TOL+1;
  d=size(Mesh.Coordinates,1);
  dmax=100000;
  while(d < dmax)
   
    iter2 = iter2+1;
    nDofs2(iter2) = size(Mesh.Coordinates,1);  
    d=nDofs2(iter2)  ;
    % Compute FE solution
  
    A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
    L = assemLoad_LFE(Mesh,P7O6(),F_HANDLE);
    [U,FreeDofs] = assemDir_LFE(Mesh,-1,GD_HANDLE);
    L = L - A*U;
    if(~isempty(FreeDofs))
      U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    end
    
    % Real and estimated error
    QuadRule=Duffy(TProd(gauleg(1,0,10)));
    Err2(iter2) = H1SErr_LFE(Mesh,U,QuadRule,UEX_HANDLE);  
    % Refine mesh by largest edge bisection
    
    tmp_Mesh = Mesh;
    Mesh = refine_REG(Mesh);
    
    % Update mesh data structure
    
    Mesh = add_Edges(Mesh);
    Mesh = rmfield(Mesh,'BdFlags');
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
    Mesh.BdFlags(Loc) = -1;   
  
    % Print out information
    
    fprintf('Eta  :  %2.4e,  nDofs  :  %2.4e\n',Err2(iter2),d);
    
    % Print out mesh if not to big
    if (nDofs2<750)
     figure;
     plot_Mesh(Mesh,'as');
     file = ['Mesh_reg_' int2str(nDofs2(iter2)) '.eps'];
     print('-depsc',file);
    end
  end
  
  save data2 Err2 iter2 nDofs2
  % Generate figures
    
  fig = figure('Name','Recovery based error estimator');
  plot(nDofs,Err,'b-x',nDofs1,Err1,'r-*',nDofs2,Err2,'g-^');
  grid('on');
  set(gca,'XScale','log','YScale','log');
  xlabel('{\bf Dofs }');
  ylabel('{\bf Errors}');
  legend('xxxxxxxxxxxxxxx','yyyyyyyyyyyyyyy','zzzzzzzzzzzzzzz','Location','SouthWest');
  p = polyfit(log(nDofs((iter-10):iter)),log(Err((iter-10):iter)),1);
  add_Slope(gca,'NorthEast',p(1));
  p = polyfit(log(nDofs1((iter1-10):iter1)),log(Err1((iter1-10):iter1)),1);
  add_Slope(gca,'SouthEast',p(1));
  p = polyfit(log(nDofs2((iter2-10):iter2)),log(Err1((iter2-10):iter2)),1);
  add_Slope(gca,'East',p(1));
   print('-depsc','adap.eps');
  % Clear memory
  
  clear all;
  
  