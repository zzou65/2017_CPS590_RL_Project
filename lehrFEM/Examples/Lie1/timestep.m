% Run script for the theta scheme (low resolution implementation). 

% Copyright 2006-2006 Patrick Meury & Mengyu Wang
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize constants
  
  NREFS = 4;        % Number of unifrom red refinements
  NSTEPS = 20;     % Number of time steps
  NFRAMES = 100;    % Number of frames
  NROUNDS = 2;      % Number of rounds
  T = 1;            % Final time
  THETA = 1;      % Coefficient of theta scheme
  G_HANDLE = @(x,t,varargin)cos(2*pi*t)*[1-x(:,1) 1-x(:,2)];    % Dirichlet boundary data
  F_HANDLE = @(x,t,varargin)-2*pi*sin(2*pi*t)*[1-x(:,1) 1-x(:,2)];    % Right hand side source term        
  U0_HANDLE = @(x,varargin)[1-x(:,1) 1-x(:,2)];;  % Initial data
  Mu_Handle=@(x,varargin)1;
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat'); 
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);  
  
  % Open up FIFO buffer
 
  buf = open(buffer(),'/scratch/users/hheumann');
    
  % Do NREFS uniform refinement steps
    
  for i = 1:NREFS
    Mesh = refine_REG(Mesh);
  end

  % Precompute matrices
  
  [I1,J1,M] = assemMat_W1F(Mesh,@MASS_W1F,Mu_Handle,P7O6());
  [I2,J2,A] = assemMat_W1F(Mesh,@STIMA_Curl_W1F,Mu_Handle,P7O6());
  
  dt = T/NSTEPS;
  S1 = sparse([I1; I2],[J1; J2],[M; dt*THETA*A]);
  S2 = sparse([I1; I2],[J1; J2],[M; -dt*(1-THETA)*A]);
  M = sparse(I1,J1,M);
  
  % Compute initial data
  
  u_old = assemLoad_W1F(Mesh,P7O6(),U0_HANDLE);
  u_old = M\u_old;
  buf = push(buf,u_old);
  
  % Start theta scheme
 
  L_old = assemLoad_W1F(Mesh,P7O6(),F_HANDLE,0);
  
  per = 0;
  progress_bar(per);
  for i = 1:NSTEPS
    if(per < floor(100*i/NSTEPS))
      per = floor(100*i/NSTEPS);
      progress_bar(per);  
    end
    
    % Assemble load vector 
    
    L_new = assemLoad_W1F(Mesh,P7O6(),F_HANDLE,i*dt);
    
    % Incorporate Dirichlet boundary data  
      
    [u_new,FreeDofs] = assemDir_W1F(Mesh,-1,G_HANDLE,gauleg(0,1,4),i*dt);
    
    % Solve the linear system
    
    rhs = S2*u_old + dt*THETA*L_new+dt*(1-THETA)*L_old;
    rhs = rhs - S1*u_new;
    u_new(FreeDofs) = S1(FreeDofs,FreeDofs)\rhs(FreeDofs);
    u_plot=plot_Norm_W1F_mov(u_new,Mesh);
    plot_Norm_W1F(u_new,Mesh);
    colorbar;
    buf = push(buf,u_plot);
    
    % Update vectors
    
    u_old = u_new;
    L_old = L_new;
    
  end
  
  % Generate movie
  
  f=movie_LFE(buf,Mesh,NFRAMES,NROUNDS);
  close(buf);
  save movie f;
  
  % Clear Memory
  
  clear all;