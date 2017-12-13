function [U,Mesh,flag,iter,ndof,err,timer,dofslvl,A,L,FreeDofs] = ...
  locmg_solve(U,Mesh,f_handle,gd_handle,errest,theta,itlvl,cyc,m,smoother,tol,maxit,varargin)
%LOCMG_SOLVE local (adaptive) multigrid solver
%   
%   [U,MESH] = LOCMG_SOLVE(U0,CMESH,F_HANDLE,GD_HANDLE,ERREST,THETA,ITLVL,CYC,M,SMOOTHER,TOL,MAXIT)
%   computes the adaptive finite element solution U of the Poisson equation
%   with right hand side given by the function handle F_HANDLE and
%   Dirichlet data given by the function handle GD_HANDLE using the local
%   (adaptive) multigrid method with longest-edge-bisection refiniments and
%   error estimator ERREST (see below) starting from the initial guess U0
%   on the (coarse) mesh CMESH.  
%
%   LOCMG_SOLVE(U0,MESH,L,A,U_BD,FDOFS,ERREST,CYC,M,SMOOTHER,TOL,MAXIT)
%   computes the finite element solution U of the equation A*U=L without
%   refining the mesh.  The mesh MESH must contain the same multilevel
%   infirmation as the returned mesh, see below.  U_BD is a vector
%   containing the boundary data for the solution U.  FDOFS specifies the
%   non-Dirichlet-boundary degrees of freedom.  ERREST is a function handle
%   taking the arguments U1,U0,A,L that returns an approximation of the
%   truncation error, where U0 and U1 are the approximate solutions before
%   and after a V-cycle.  This is used as a stopping criterium.  If ERREST 
%   is empty, the energy norm of U1-U0 is used; if ERREST is a vector U_E, 
%   then the energy norm of U1-U_E is used.  When using this syntax, some
%   of the descriptions below no longer apply.
%
%   The initial guess U0 may be replaced by an empty matrix, in which case
%   the zero vector of the correct size is used.
%
%   M denotes the number of smoothing steps.  If M is a scalar, then M
%   presmoothing steps and M postsmoothing steps are used.  If M=[M1,M2],
%   then the algorithm uses M1 presmoothing steps and M2 postsmoothing
%   steps.
%
%   CYC denotes the number of coarse grid corrections to be used at each
%   multigrid iteration, ie. CYC=1 for V-cylces and CYC=2 for W-cycles.
%
%   SMOOTHER is a function handle for the smoother.  It must take the
%   arguments x,A,b in that order.  To pass further arguments to the
%   smoother, use LOCMG_SOLVE(...,VARARGIN).
%
%   ERREST determines the aposteriori error estimetor.  It can take any
%   one of the following values:
%    'r' residual-based error estimator.
%    'c' recovery-based error estimator.
%    'h' hierarchical error estimator.
%   Alternatively, ERREST can be a function handle of the form
%       ETA = ERREST(U,MESH,A)
%   where U is the approximate solution, MESH is the current mesh and ETA
%   is the estimated error.  ETA is a vector of length M if the Mesh has M
%   elements which contains the error on each element of the mesh.  The
%   elements with ETA >= THETA*max(ETA) are refined.
%
%   THETA determines the sensitivity of the refinements.  After every
%   iteration, cells with error (estimated by ERREST) greater than THETA
%   times the maximal error are refined.
%
%   ITLVL is the number of multigrid cycles that are done in each
%   iteration, i.e. between refinements.
%
%   The tolerance TOL and the maximum number of iterations MAXIT are used
%   as stopping criteria for the multigrid iterations.  Reasonable values
%   of these arguments are TOL = 0.05 and MAXIT = 30.
%
%   If the argument F_HANDLE is a function handle packed in a cell array,
%   then it is treated as the full assembly routine for the load vector.
%   It must take the single argument Mesh.
%
%   The returned mesh MESH contains the following fields (in addition to
%   the usual fields) : (N : number of vertices, M : number of elements)
%     VertLevel   Nx1 array containig the lowest refinement level where a
%                 vertex appears (0 is the coarsest level)
%     ElemLevel   Mx1 array containing the lowest refinement level where 
%                 all vertices of an element are in the mesh
%     VertInLvl   Nx1 array containing the highest refinement level for
%                 wich a vertex is inside the refinement zone, ie the
%                 minimal level of an adjacent element
%     VertClLvl   Nx1 array containing the highest refinement level for
%                 which a vertex is in the closure of the refinement zone,
%                 ie the maximal level of an adjacent element
%     FatherVert  Nx2 array containg the indices of the endpoints of the
%                 edge on which a vertex was created through longest edge
%                 bisection
%   If the input mesh CMESH already contains these fields, then the initial
%   multigrid structure is defined accordingly; otherwise, all the elements
%   of CMESH are considered to be on the coarsest level 0.
%
%   [U,MESH,FLAG] = LOCMG_SOLVE(...) returns the logical FLAG, which
%   specifies whether or not the method converged within MAXIT steps to 
%   the prescribed tolerance TOL according to the error estimator ERREST.
%
%   [U,MESH,FLAG,ITER,NDOF,ERR]  = LOCMG_SOLVE(...) also returns
%    ITER  the number of calculated iterations.
%    NDOF  a 1 by ITER array containg the number of degrees of feedom at
%          each iteration.
%    ERR  a 1 by ITER array containing the error estimate at each level.
%
%   [U,MESH,FLAG,ITER,NDOF,ERR,TIMER]  = LOCMG_SOLVE(...) also returns the
%   time in seconds taken up by various parts of the code.  TIMER is a
%   struct with the following fields:
%    total   total time
%    mg      time required for multigrid cycles
%    errest  time required for error estimator 
%    mesh    time required for mesh generation
%    stima   time required for assembly of stiffness matrix and load vector
%    ref     time required for mesh refinement
%
%   [U,MESH,FLAG,ITER,NDOF,ERR,TIMER,DOFSLVL]  = LOCMG_SOLVE(...) also
%   returns the LVLxITER matrix DOFSLVL, which contains the number of
%   degrees of freedom on each level at each iteration.  DOFSLVL(lvl+1,i) is
%   the number of degrees of freedom smoothed on level lvl at iteration i.
%
%   [U,MESH,FLAG,ITER,NDOF,ERR,TIMER,DOFSLVL,A,L,FDOFS]  = LOCMG_SOLVE(...) also
%   returns the stiffness matrix A on the full mesh, the corresponding
%   load vector L amd the free degrees of freedom FDOFS.
%
%   Example:
%
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%   U0 = zeros(size(Mesh.Coordinates,1),1);
%   [U,Mesh,flag,iter,ndof,err,timer] = ...
%   locmg_solve(U0,Mesh,@f_LShap,@g_D_LShap,'c',0.5,1,1,[2,1],@gs_smooth,0.05,20);
%
%   See also multigrid_solve, ErrEst_RES, ErrEst_REC, ErrEst_HIER,
%   refine_LEB, assemLoad_LFE_S1.

%   Copyright 2006-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland


  % Check that U has correct size
  
  if(~all(size(U)==[size(Mesh.Coordinates,1),1]))
    if(~isempty(U))
      fprintf('The initial guess has dimension %dx%d instead of %dx%d.\n  Replacing U by zero vector...\n',size(U,1),size(U,2),size(Mesh.Coordinates,1),1);
    end
    U = zeros(size(Mesh.Coordinates,1),1);
  end

  % Initialize timers
  
  timer.total = 0;
  timer.mg = 0;
  timer.errest = 0;
  timer.mesh = 0;
  timer.stima = 0;
  timer.ref = 0;
  
  t_tot = cputime;
  
  % Initialize constants

  QUADRULE = P7O6();
  
  % Define assembly routine for load vector
  
  refine_mesh = true;
  if(isa(f_handle,'cell'))
    assemLoad = f_handle{1};
  elseif(isa(f_handle,'numeric')) % alternative syntax without mesh refinements
    refine_mesh = false;
  else
    assemLoad = @(mesh) assemLoad_LFE(mesh,QUADRULE,f_handle);
  end

  % Initialize mesh, assuming that if it contains the field 'VertLevel', it
  % must contain all the other relevent fields
  
  if(~isfield(Mesh,'VertLevel'))
    t = cputime;
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);  
    Mesh = add_Edges(Mesh);                                
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
    Mesh.BdFlags(Loc) = -1;
    Mesh = init_LEB(Mesh);
    Mesh.VertLevel = zeros(size(Mesh.Coordinates,1),1);
    Mesh.FatherVert = zeros(size(Mesh.Coordinates,1),2);
    Mesh.ElemLevel = zeros(size(Mesh.Elements,1),1);
    Mesh.VertInLvl = zeros(size(Mesh.Coordinates,1),1);
    Mesh.VertClLvl = zeros(size(Mesh.Coordinates,1),1);
    Mesh = add_Patches(Mesh);
    timer.mesh = cputime-t+timer.mesh;
  end
  
  % Initialization with mesh refinements
  
  if(refine_mesh)  
  
    % Assemble stiffness matrix and load vector

    t = cputime;
    A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
    L = assemLoad(Mesh);
    [U_bd,FreeDofs] = assemDir_LFE(Mesh,-1,gd_handle);  
    
    timer.stima = cputime-t+timer.stima;

    % Check for custom error estimator

    if(isa(errest,'function_handle'))
      errest_handle = errest;
      errest = 'function_handle';
    end
    
  % Initialization without mesh refinements
    
  else
    
    % reassign arguments
    
    L = f_handle;
    A = gd_handle;
    U_bd = errest;
    FreeDofs = theta;
    errest = itlvl;
    
    itlvl = 1;
    theta = 1;
    
    % initialize error estimator
    
    if(isempty(errest))
      errest = @(u1,u0,varargin) sqrt((u1-u0)'*A*(u1-u0));
    elseif(isa(errest,'numeric'))
      u_ex = errest;
      errest = @(u,varargin) sqrt((u-u_ex)'*A*(u-u_ex));
    end
    
  end
  
  % Define free dofs
  
  FDofs = false(size(Mesh.Coordinates,1),1);
  FDofs(FreeDofs) = true;
  U(~FDofs) = U_bd(~FDofs);
  
  % Construct stiffness matrix on all levels
  
  t = cputime;
  [fA,refZone,cAind] = getA(A,Mesh.VertLevel,Mesh.VertClLvl,Mesh.FatherVert,max(Mesh.VertLevel));
  timer.stima = cputime-t+timer.stima;
  
  % Iterate
  
  flag = false;
  ndof = zeros(1,maxit);
  err = zeros(1,maxit);
  dofslevel = cell(1,maxit);

  for iter=1:maxit
    
    % Count degrees of freedom : total and on each level
    
    ndof(iter) = nnz(FDofs);
    
    LVL = max(Mesh.VertLevel);
    dofslevel{iter} = zeros(1,LVL+1);
    for lvl=0:LVL
      dofslevel{iter}(lvl+1) = nnz(Mesh.VertInLvl(FDofs)>=lvl & Mesh.VertLevel(FDofs)<=lvl);
    end
   
    % Compute solution
    %
    % if(any(FDofs))
    %   L = L - A*U_bd;
    %   U(FDofs) = A(FDofs,FDofs)\L(FDofs);
    % end
    
    t = cputime;
    for i=1:itlvl
      U0 = U;
      U = locmg_step(U0,fA,cAind,L,Mesh.VertLevel,Mesh.VertInLvl,Mesh.FatherVert,FDofs,refZone,LVL,smoother,m,cyc,varargin{:});
    end
    timer.mg = cputime-t+timer.mg;
    
    % Estimate error
    
    t = cputime;
    if(refine_mesh)
      switch errest
        case 'r'    % residual-based error estimator
          Mesh = add_Edge2Elem(Mesh);
          Eta_K = ErrEst_RES(U,Mesh,QUADRULE,f_handle);
        case 'c'    % recovery-based error estimator
          Eta_K = ErrEst_REC(U,Mesh,QUADRULE);
        case 'h'    % hierarchical error estimator
          Mesh = add_Edge2Elem(Mesh);
          R_Handle = @Res_Lapl;
          S_Handle = @STIMA_Lapl_ErrEst;
          Eta_K = ErrEst_HIER(U,Mesh,R_Handle,S_Handle,QUADRULE,f_handle);
        case 'function_handle'  % custom error estimator
          Eta_K = errest_handle(U,Mesh);
        otherwise
          error('Invalid error estimator.');
      end
      err(iter) = sqrt(sum(Eta_K.^2));
    else
      err(iter) = errest(U,U0,A,L);
    end
    timer.errest = cputime-t+timer.errest;
    
    % Exit if tolerance is reached
    
    if(err(iter) < tol)
      flag = true;
      break;
    elseif(iter < maxit && theta <= 1 && refine_mesh)
      
      t = cputime;
    
      % Mark elements for refinement
    
      Eta_max = max(Eta_K);
      MarkedElem = find(Eta_K >= Eta_max*theta);

      % Refine mesh by largest edge bisection
      
      cNumVert = size(Mesh.Coordinates,1);
      
      [Mesh,fatherVert] = refine_LEB(Mesh,MarkedElem);
      
      timer.ref = cputime-t+timer.ref;
      
      
      t = cputime;

      % Update mesh data structure

      Mesh = add_Edges(Mesh);
      Mesh = add_Patches(Mesh);

      numVert = size(Mesh.Coordinates,1);
      numElem = size(Mesh.Elements,1);     

      % Add multilevel information to mesh :
      %   VertLevel : the lowest refinement level where the vertex appears
      %   ElemLevel : the lowest refinement level where all vertices of an
      %               element are in the mesh
      %   VertInLvl : the highest refinement level for which the vertex is
      %               inside of the refinement zone
      %   VertClLvl : the highest refinement level for which the vertex is
      %               in the closure of the refinement zone.

      Mesh.FatherVert = [Mesh.FatherVert;fatherVert];

      Level = Mesh.VertLevel;
      Mesh.VertLevel = zeros(numVert,1);
      Mesh.VertLevel(1:size(Level,1)) = Level;
      for i=cNumVert+1:numVert
        Mesh.VertLevel(i) = max(Mesh.VertLevel(Mesh.FatherVert(i,:)))+1;
      end

      Mesh.ElemLevel = zeros(numElem,1);
      for i=1:numVert
        Mesh.ElemLevel(Mesh.AdjElements(i,1:Mesh.nAdjElements(i))) = ...
          max(Mesh.ElemLevel(Mesh.AdjElements(i,1:Mesh.nAdjElements(i))),...
              Mesh.VertLevel(i));
      end

      Mesh.VertInLvl = zeros(numVert,1);
      Mesh.VertClLvl = zeros(numVert,1);
      for i=1:numVert
        Mesh.VertInLvl(i) = min(Mesh.ElemLevel(Mesh.AdjElements(i,1:Mesh.nAdjElements(i))));
        Mesh.VertClLvl(i) = max(Mesh.ElemLevel(Mesh.AdjElements(i,1:Mesh.nAdjElements(i))));
      end

      % Prolong solution to new mesh

      U = [U;zeros(numVert-cNumVert,1)];
      for i=cNumVert+1:numVert
        U(i) = 0.5*(sum(U(Mesh.FatherVert(i,:))));
      end
%       U(cNumVert+1:numVert) = 0.5*sum(U(Mesh.FatherVert(cNumVert+1:numVert,:)).').';


      % Find boundary edges

      Mesh = rmfield(Mesh,'BdFlags');
      Loc = get_BdEdges(Mesh);
      Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
      Mesh.BdFlags(Loc) = -1;

      % Handle Dirichlet boundary conditions

      [U_bd,FreeDofs] = assemDir_LFE(Mesh,-1,gd_handle);  
      FDofs = false(size(Mesh.Coordinates,1),1);
      FDofs(FreeDofs) = true;
      U(~FDofs) = U_bd(~FDofs);

      timer.mesh = cputime-t+timer.mesh;

      t = cputime;

      % Construct mesh in refinement zone

      newVert = [false(1,cNumVert),true(1,numVert-cNumVert)]';
      refVert = false(1,numVert)';
      refElem = false(1,numElem)';
      for i=cNumVert+1:numVert
        refElem(Mesh.AdjElements(i,1:Mesh.nAdjElements(i))) = true;
      end
      refVert(Mesh.Elements(refElem,:)) = true;
      oldVert = refVert & ~newVert;

      rNewVert = newVert(refVert);
      rOldVert = oldVert(refVert);

      f2r = cumsum(refVert);  % full index to refinement area index
      f2r(~refVert) = 0;

      if(exist('rMesh','var'))
        clear('rMesh');
      end

      rMesh.Coordinates = Mesh.Coordinates(refVert,:);
      rMesh.Elements = f2r(Mesh.Elements(refElem,:));
      rMesh.ElemFlag = ones(size(rMesh.Elements,1),1);

      timer.mesh = cputime-t+timer.mesh;

      % Assemble stiffness matrix and load vector in refinement zone

      t = cputime;

      rA = assemMat_LFE(rMesh,@STIMA_Lapl_LFE);
      rL = assemLoad(rMesh);

      % Construct full load vector

      L = [L;zeros(numVert-cNumVert,1)];
      L(newVert) = rL(rNewVert);
      ind = Mesh.FatherVert(newVert,:);
      L(ind(:,1)) = L(ind(:,1)) - 0.5*rL(rNewVert);
      L(ind(:,2)) = L(ind(:,2)) - 0.5*rL(rNewVert);

      %% Construct full stiffness matrix
      
      % update entries of A

      A(numVert,numVert) = 0; % increase size of A

      A(newVert,refVert) = rA(rNewVert,:);
      A(oldVert,newVert) = rA(rOldVert,rNewVert);
      
      % extract relevant entries of A
      
      [A_I,A_J,A_A] = find(A);

      t_ind = newVert(A_I);
      I = A_I(t_ind);
      J = A_J(t_ind);
      ta = A_A(t_ind);
      
      % prealocate memory for correction
      
      len = length(A_A);
      d_len = 4*length(I);
      A_I = [A_I;zeros(d_len,1)];
      A_J = [A_J;zeros(d_len,1)];
      A_A = [A_A;zeros(d_len,1)];
      
      % calculate entries of correction
      
      for k=1:length(I)
        if(oldVert(J(k)))
          i = Mesh.FatherVert(I(k),:);
          A_I(len+1:len+2) = i;
          A_J(len+1:len+2) = J(k);
          A_I(len+3:len+4) = J(k);
          A_J(len+3:len+4) = i;
          A_A(len+1:len+4) = -0.5*ta(k);
          len = len + 4;
        else%if(newVert(J(k)))
          i = Mesh.FatherVert(I(k),:);
          j = Mesh.FatherVert(J(k),:);
          A_I(len+1:len+2) = i;
          A_I(len+3:len+4) = i;
          A_J(len+1:2:len+3) = j;
          A_J(len+2:2:len+4) = j;
          A_A(len+1:len+4) = -0.25*ta(k);
          len = len + 4;
        end
      end
      
      % assemble full stiffness matrix
      
      A = sparse(A_I,A_J,A_A,numVert,numVert);

      timer.stima = cputime-t+timer.stima; 

      % Construct stiffness matrix on all levels
      
      t = cputime;
    
      [fA,refZone,cAind] = getA(A,Mesh.VertLevel,Mesh.VertClLvl,Mesh.FatherVert,max(Mesh.VertLevel));

      timer.stima = cputime-t+timer.stima;

    end
  end

  % Define return arguments
  
  ndof = ndof(1:iter);
  err = err(1:iter);
%   dofslevel = {dofslevel{1:iter}};
  dofslvl = zeros(LVL+1,iter);
  for i=1:iter
    dofslvl(1:length(dofslevel{i}),i) = dofslevel{i};
  end
  
  timer.total = cputime-t_tot+timer.total;
  
return



%%% Construct stiffness matrix on all levels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function [fA,refZone,cAind] = getA(A,VertLevel,VertClLvl,FatherVert,LVL)
% GETA construct and store stiffness matrix on coarse levels
%
%   [FA,REFZONE,CAIND] = GETA(A,VERTLEVEL,VERTCLLVL,FATHERVERT,LVL)
%   constructs the stiffness matrix and relevent indices on all levels from
%   the stiffness matrix A on the full mesh.
%
%   VERTLEVEL is an N = size(A,1) = size(A,2) by 1 vector containing the
%   level (0 to LVL) of each vertex.  The level is the lowest refinement
%   level in which the vertex appears.  VERTCLLVL contains the highest
%   level for which the vertex is in the closure of the refinement uone, ie
%   the maximal level of the adjacent elements.
%
%   FATHERVERT is an N by 2 matrix.  For every vertex i that has positive
%   level, ie is not in the coarsest grid, FATHERVERT(i,:) are the
%   endpoints of the edge on which the vertex i was constructed through
%   longest edge bisection.
%
%   The return argument FA is a 1 by LVL+1 cell array containing the
%   necessary elements of the stiffness matrix on each level, ie, on level
%   lvl, the rows corresponding to vertices with VERTCLLVL >= lvl.  These
%   are given by the logical index REFZONE{lvl}.
%
%   CAIND{lvl} is a logical index of size size(FA{lvl},1) which restricts
%   the stiffness matrix to the refinement zone of the next-finer level.
%
%   See also locmg_step, locmg_solve.

%   Copyright 2006-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  fA = cell(1,LVL+1);
  refZone = cell(1,LVL);
  crefZone = cell(1,LVL);
  cAind = cell(1,LVL);

  for lvl=LVL:-1:1

    % Determine fine and coarse nodes

    fVertFlag = VertLevel == lvl;
    cVertFlag = ~fVertFlag;

    refZone{lvl} = VertClLvl >= lvl;
    
    crefZone{lvl} = refZone{lvl}(cVertFlag);
    if(lvl < LVL)
      cAind{lvl+1} = crefZone{lvl+1}(refZone{lvl});
    end

    f2c = cumsum(cVertFlag);  % lvl-index to lvl-1-index, 0 on fine grid points
    f2c(fVertFlag) = 0;

    % Store stiffness matrix in refinement zone

    fA{lvl+1} = A(refZone{lvl},:);   

    %% Construct coarse-grid stiffness matrix
    
    % extract relevant entries of A
    
    [A_I,A_J,A_A] = find(A);
    ci = cVertFlag(A_I);
    cj = cVertFlag(A_J);
    ind = ci & cj;
    
    t_ind = ~ci;
    I = A_I(t_ind);
    J = A_J(t_ind);
    ta = A_A(t_ind);
    
    % prealocate memory for correction

    maxlen = 4*length(I);
    dA_I = zeros(maxlen,1);
    dA_J = zeros(maxlen,1);
    dA_A = zeros(maxlen,1);
    len = 0;
    
    % calculate entries of correction
    
    for k=1:length(I)
      if(cVertFlag(J(k)))
        i = FatherVert(I(k),:);
        dA_I(len+1:len+2) = i;
        dA_J(len+1:len+2) = J(k);
        dA_I(len+3:len+4) = J(k);
        dA_J(len+3:len+4) = i;
        dA_A(len+1:len+4) = 0.5*ta(k);
        len = len + 4;
      else
        i = FatherVert(I(k),:);
        j = FatherVert(J(k),:);
        dA_I(len+1:len+2) = i;
        dA_I(len+3:len+4) = i;
        dA_J(len+1:2:len+3) = j;
        dA_J(len+2:2:len+4) = j;
        dA_A(len+1:len+4) = 0.25*ta(k);
        len = len + 4;
      end
    end
    
    % assemble coarse-grid stiffness matrix
    
    cdi = cVertFlag(dA_I);
    cdj = cVertFlag(dA_J);
    d_ind = cdi & cdj;
    
    I = f2c([A_I(ind);dA_I(d_ind)]);
    J = f2c([A_J(ind);dA_J(d_ind)]);
    a = [A_A(ind);dA_A(d_ind)];
    
    A = sparse(I,J,a);
    
    % Redefine variables for next level

    VertLevel = VertLevel(cVertFlag);
    VertClLvl = VertClLvl(cVertFlag);
    FatherVert = FatherVert(cVertFlag,:);
    FatherVert(FatherVert~=0) = f2c(FatherVert(FatherVert~=0));

  end

  % Store stiffness matrix in coarsest mesh
  
  fA{1} = A;
  
  if(LVL>0)
    cAind{1} = crefZone{1};
  end

return
  
  
  
%%% Multigrid iteration %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
function U = locmg_step(U,fA,cAind,L,VertLevel,VertInLvl,FatherVert,FDofs,refZone,lvl,smoother,m,cyc,varargin)
% LOCMG_STEP local (adaptive) multigrid step
%
%   U = LOCMG_STEP(U,FA,CAIND,L,VERTLEVEL,VERTINLVL,FATHERVERT,FDOFS,REFZONE,LVL,SMOOTHER,M,CYC)
%   performs a multigrid cycle starting with approximation U, with
%   stiffness matrix A on the finest level and load vector L on the finest
%   level.  Note that Dirichlet data should not yet be integrated into L
%   and the input arguments should have dimension equal to the
%   full number N of vertices.
%
%   M(1) presmoothing steps and M(end) postsmoothing steps are used.  CYC
%   denotes the number of coarse grid corrections, ie CYC=1 for V-cycles
%   and CYC=2 for W-cycles.
%
%   SMOOTHER is a function handle for the smoother.  It must take the
%   arguments x,A,b in that order.  To pass further arguments to the
%   smoother, use LOCMG_STEP(...,VARARGIN).
%
%   The stiffness matrix is given by the cell array FA.  More precicely,
%   FA{lvl+1} is the REFZONE{lvl} rows of the stiffness matrix on level lvl,
%   where REFZONE{lvl} is a logical index of the appropriate dimension.
%   CAIND{lvl} is a logical index of size size(FA{lvl},1) which restricts
%   the stiffness matrix to the refinement zone of the next-finer level.
%   These three cell arrays can be constructed with the getA routine.
%
%   VERTLEVEL is an N by 1 vector containing the level (0 to LVL) of each
%   vertex.  The level is the lowest refinement level in which the vertex
%   appears.  VERTINLVL contains the highest level for which the vertex is
%   inside the refinement zone, ie the minimal level of the adjacent
%   elements.  Degrees of freedom with VERTINLVL >= LVL and which are not
%   on the boundary are smoothed.
%
%   FATHERVERT is an N by 2 matrix.  For every vertex i that has positive
%   level, ie is not in the coarsest grid, FATHERVERT(i,:) are the
%   endpoints of the edge on which the vertex i was constructed through
%   longest edge bisection.
%
%   FDOFS is an N by 1 logical array with true at the free degrees of
%   freediom and false at the boundary vertices.
%
%   See also getA, locmg_solve.

%   Copyright 2006-2006 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland


  % Use direct solver on the coarsest mesh

  if(lvl == 0)
    if(any(FDofs))
      L(FDofs) = L(FDofs) - fA{1}(FDofs,~FDofs)*U(~FDofs);
      U(FDofs) = fA{1}(FDofs,FDofs)\L(FDofs);
    end
    return
  end

  % Determine degrees of freedom and fine and coarse nodes
  
  fVertFlag = VertLevel == lvl;
  cVertFlag = ~fVertFlag;
  
  fVertFlagRef = fVertFlag(refZone{lvl});
  cVertFlagRef = cVertFlag(refZone{lvl});
  
  dofs = (VertInLvl >= lvl & FDofs); 
  dofsRef = dofs(refZone{lvl});
  notdofs = ~dofs;
  
  cRefZone = refZone{lvl}(cVertFlag);
  
  % Relax M(1) times 

  if(any(dofs) && m(1)>0)
    L_dofs = L(dofs)-fA{lvl+1}(dofsRef,notdofs)*U(notdofs);
    for i = 1:m(1)
      U(dofs) = smoother(U(dofs),fA{lvl+1}(dofsRef,dofs),L_dofs,varargin{:});
    end
  end
  
  % Initialize coarse grid correction

  f2c = cumsum(cVertFlag);  % lvl-index to lvl-1-index, 0 on fine grid points
  f2c(fVertFlag) = 0;
  
  % Calculate and restrict residual
  
  res = L(refZone{lvl}) - fA{lvl+1}*U;                % residual
  ind = f2c(FatherVert(fVertFlag,:).').';
  cRes = zeros(size(fA{lvl},2),1);                    % restricted residual 
  cRes(cRefZone) = res(cVertFlagRef);
  cRes(ind(:,1)) = cRes(ind(:,1)) + 0.5*res(fVertFlagRef);
  cRes(ind(:,2)) = cRes(ind(:,2)) + 0.5*res(fVertFlagRef);
  cRes = cRes(cRefZone);
  
  % Calculate residual term for hU and its restriction to the coarse grid
  
  cDelta = L(cVertFlag);
  cDelta(cRefZone) = fA{lvl}(cAind{lvl},:)*U(cVertFlag) + cRes;
  
  % Initialize coarse grid correction
  
  cVertLevel = VertLevel(cVertFlag);
  cVertInLvl = VertInLvl(cVertFlag);
  cFatherVert = FatherVert(cVertFlag,:);
  cFatherVert(cFatherVert~=0) = f2c(cFatherVert(cFatherVert~=0));
  cFDofs = FDofs(cVertFlag);
  
  % Do coarse grid correction 
  
  U_old = sum(U(FatherVert(fVertFlag,:)).').';
  for i=1:cyc
    U(cVertFlag) = locmg_step(U(cVertFlag),fA,cAind,cDelta,cVertLevel,cVertInLvl,cFatherVert,cFDofs,refZone,lvl-1,smoother,m,cyc,varargin{:});
  end
  U(fVertFlag) = U(fVertFlag) + 0.5*(sum(U(FatherVert(fVertFlag,:)).').'-U_old);

  % Relax M(2) times 

  if(any(dofs) && m(end)>0)
    L_dofs = L(dofs)-fA{lvl+1}(dofsRef,notdofs)*U(notdofs);
    for i = 1:m(end)
      U(dofs) = smoother(U(dofs),fA{lvl+1}(dofsRef,dofs),L_dofs,varargin{:});
    end
  end
  
return