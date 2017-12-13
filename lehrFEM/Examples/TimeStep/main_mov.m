% Run script for the theta scheme (low resolution implementation). 

% Copyright 2006-2006 Patrick Meury & Mengyu Wang
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize constants
  
  NREFS = 6;        % Number of unifrom red refinements
  NSTEPS = 200;     % Number of time steps
  NFRAMES = 100;    % Number of frames
  NROUNDS = 2;      % Number of rounds
  T = 1;            % Final time
  THETA = 0.5;      % Coefficient of theta scheme
  G_HANDLE = @g;    % Dirichlet boundary data
  F_HANDLE = @f;    % Right hand side source term        
  U0_HANDLE = @u0;  % Initial data
  
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat'); 
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);  
  
  % Open up FIFO buffer
 
  buf = open(buffer(),'/scratch/users/meury');
    
  % Do NREFS uniform refinement steps
    
  for i = 1:NREFS
    Mesh = refine_REG(Mesh);
  end

  % Precompute matrices
  
  [I1,J1,M] = assemMat_LFE(Mesh,@MASS_LFE);
  [I2,J2,A] = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
  
  dt = T/NSTEPS;
  S1 = sparse([I1; I2],[J1; J2],[M; dt*THETA*A]);
  S2 = sparse([I1; I2],[J1; J2],[M; -dt*(1-THETA)*A]);
  M = sparse(I1,J1,M);
  
  % Compute initial data
  
  u_old = assemLoad_LFE(Mesh,P3O3(),U0_HANDLE);
  u_old = M\u_old;
  buf = push(buf,u_old);
  
  % Start theta scheme
 
  L_old = assemLoad_LFE(Mesh,P3O3(),F_HANDLE,0);
  
  per = 0;
  progress_bar(per);
  for i = 1:NSTEPS
    if(per < floor(100*i/NSTEPS))
      per = floor(100*i/NSTEPS);
      progress_bar(per);  
    end
    
    % Assemble load vector 
    
    L_new = assemLoad_LFE(Mesh,P3O3(),F_HANDLE,i*dt);
    
    % Incorporate Dirichlet boundary data  
      
    [u_new,FreeDofs] = assemDir_LFE(Mesh,-1,G_HANDLE,i*dt);
    
    % Solve the linear system
    
    rhs = S2*u_old + dt*THETA*L_new+dt*(1-THETA)*L_old;
    rhs = rhs - S1*u_new;
    u_new(FreeDofs) = S1(FreeDofs,FreeDofs)\rhs(FreeDofs);
    buf = push(buf,u_new);
    
    % Update vectors
    
    u_old = u_new;
    L_old = L_new;
    
  end
  
  % Generate movie
  
  movie_LFE(buf,Mesh,NFRAMES,NROUNDS);
  close(buf);
  
  % Clear Memory
  
  clear all;
  