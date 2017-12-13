% Run script for the theta scheme (high resolution implementation). 

% Copyright 2006-2006 Patrick Meury & Mengyu Wang
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize constants

  IREFS = 2;        % Number of initial mesh refinements [Do not change]
  NREFS = 6;        % Number of unifrom red refinements
  NSTEPS = 200;     % Number of time steps
  T = 1;            % Final time
  THETA = 0.5;      % Coefficient of theta scheme
  G_HANDLE = @g;    % Dirichlet boundary data
  F_HANDLE = @f;    % Right hand side source term        
  U0_HANDLE = @u0;  % Initial data
  TOL = 1e-12;      % Stopping criterion
  MAXIT = 100;      % Maximum number of iterations
  
  % Initialize mesh
  
  CMesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  CMesh = add_Edges(CMesh);
  Loc = get_BdEdges(CMesh);
  CMesh.BdFlags = zeros(size(CMesh.Edges,1),1);
  CMesh.BdFlags(Loc) = -1;
  CMesh.ElemFlag = zeros(size(CMesh.Elements,1),1);
  for i = 1:IREFS
    CMesh = refine_REG(CMesh);
  end
  Loc = get_BdEdges(CMesh);
  DDofs = unique([CMesh.Edges(Loc,1); CMesh.Edges(Loc,2)]);
  CFDofs = setdiff(1:size(CMesh.Coordinates,1),DDofs);
  
  % Compute initial data and generate multilevel data structures
  
  dt = T/NSTEPS;
  [I1,J1,A] = assemMat_LFE(CMesh,@STIMA_Lapl_LFE);
  [I2,J2,M] = assemMat_LFE(CMesh,@MASS_LFE);
  L_old = assemLoad_LFE(CMesh,P1O2(),U0_HANDLE);
  U_old = zeros(size(L_old));  
  for i = 1:NREFS       
    S1 = sparse([I1; I2],[J1; J2],[M; dt*THETA*A]);
    BPX_STIMA.D{i} = diag(S1(CFDofs,CFDofs));
    M = sparse(I2,J2,M);
    BPX_MASS.D{i} = diag(M);
    
    % Compute solution on the current level
    
    U_old = pcg_solve(U_old,M,L_old,TOL,MAXIT,@bpx_prec,BPX_MASS);
    
    % Refine the mesh
    
    FMesh = refine_REG(CMesh);
    Loc = get_BdEdges(FMesh);
    DDofs = unique([FMesh.Edges(Loc,1); FMesh.Edges(Loc,2)]);
    FFDofs = setdiff(1:size(FMesh.Coordinates,1),DDofs);
    
    % Compute stiffness matrix and load vector
    
    [I1,J1,A] = assemMat_LFE(FMesh,@STIMA_Lapl_LFE);
    [I2,J2,M] = assemMat_LFE(FMesh,@MASS_LFE);
    L_old = assemLoad_LFE(FMesh,P1O2(),U0_HANDLE); 
    
    % Compute prolongation matrix
    
    P = get_PMat_LFE(CMesh,FMesh);
    BPX_MASS.P{i} = P;
    BPX_STIMA.P{i} = P(FFDofs,CFDofs);
    
    % Update coarse mesh
    
    CMesh = FMesh;
    CFDofs = FFDofs;
    
    % Prolongate initial value
    
    U_old = P*U_old;
    
  end  
  clear('FMesh','P');
  
  % Compute load vector and stiffness matrix on finest mesh
  
  [I1,J1,A] = assemMat_LFE(CMesh,@STIMA_Lapl_LFE);
  [I2,J2,M] = assemMat_LFE(CMesh,@MASS_LFE);
  L_old = assemLoad_LFE(CMesh,P1O2(),U0_HANDLE);
  S1 = sparse([I1; I2],[J1; J2],[M; dt*THETA*A]);
  S2 = sparse([I1; I2],[J1; J2],[M; -dt*(1-THETA)*A]);
  BPX_STIMA.D{NREFS+1} = diag(S1(CFDofs,CFDofs));
  M = sparse(I2,J2,M);
  BPX_MASS.D{NREFS+1} = diag(M);
  clear('I1','J1','A','I2','J2')
  
  % Run preconditioned CG solver
  
  [U_old,flag] = pcg_solve(U_old,M,L_old,TOL,MAXIT,@bpx_prec,BPX_MASS);
  if(flag == 0)
    error('PCG did not converge');
  end
  clear('BPX_MASS','M');
  
  % Run theta scheme
  
  per = 0;
  progress_bar(per);
  
  L_old = assemLoad_LFE(CMesh,P1O2(),F_HANDLE,0);
  for i = 1:NSTEPS
    if(per < floor(100*i/NSTEPS))
      per = floor(100*i/NSTEPS);
      progress_bar(per);  
    end
    
    % Assemble load vector 
    
    L_new = assemLoad_LFE(CMesh,P1O2(),F_HANDLE,i*dt);
    
    % Incorporate Dirichlet boundary data  
      
    [U_new,FreeDofs] = assemDir_LFE(CMesh,-1,G_HANDLE,i*dt);
    
    % Solve the linear system using PCG with BPX preconditioner
    
    rhs = S2*U_old + dt*THETA*L_new+dt*(1-THETA)*L_old;
    rhs = rhs - S1*U_new;
    [U_new(FreeDofs),flag] = pcg_solve(U_old(FreeDofs),S1(FreeDofs,FreeDofs),rhs(FreeDofs), ...
                                       TOL,MAXIT,@bpx_prec,BPX_STIMA);
    if(flag == 0)
      error('PCG did not converge');  
    end
    
    % Update vectors
    
    U_old = U_new;
    L_old = L_new;
    
  end
  
  % Generate figure
  
  plot_LFE(U_new,CMesh);
  colorbar;
  
  % Clear memory
  
  clear all;
  