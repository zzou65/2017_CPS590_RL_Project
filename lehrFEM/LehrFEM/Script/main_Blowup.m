% Visualizes instabilities of the explicit Euler scheme when the CFL
% condition is violated. 

% Copyright 2006-2006 Patrick Meury & Mengyu Wang
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize constants
  
  NREFS = 4;                                                                                           % Number of unifrom red refinements 
  NSTEPS = 3200;                                                                                       % Number of time steps
  T = 1;                                                                                               % Final time
  G = @(x,varargin)zeros(size(x,1),1);                                                                 % Dirichlet boundary data
  F = @(x,ElemFlag,t,varargin)2*pi*(4*pi*cos(2*pi*t)-sin(2*pi*t))*sin(2*pi*x(:,1)).*sin(2*pi*x(:,2));  % Right hand side source term        
  U0 = @(x,varargin)sin(pi*x(:,1)).*sin(pi*x(:,2));                                                    % Initial data
  DUMMY = @(x,varargin)zeros(size(x,1),1);                                                             % Dummy routine
  
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat'); 
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = -1;
  Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);  
    
  % Do NREFS uniform refinement steps
    
  for i = 1:NREFS
    Mesh = refine_REG(Mesh);
  end

  % Precompute matrices
  
  [I1,J1,M] = assemMat_LFE(Mesh,@MASS_LFE);
  [I2,J2,A] = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
  
  dt = T/NSTEPS;
  S = sparse([I1; I2],[J1; J2],[M; -dt*A]);
  M = sparse(I1,J1,M);
  
  % Compute initial data
  
  u_old = assemLoad_LFE(Mesh,P1O2(),U0);
  u_old = M\u_old;
  
  L2Norm = zeros(1,NSTEPS+1);
  H1SNorm = zeros(1,NSTEPS+1);
  LInfNorm = zeros(1,NSTEPS+1);
  L2Norm(1) = L2Err_LFE(Mesh,u_old,P1O2(),DUMMY);
  H1SNorm(1) = H1SErr_LFE(Mesh,u_old,P1O2(),DUMMY);
  LInfNorm(1) = LInfErr_LFE(Mesh,u_old,DUMMY);
  
  % Start theta scheme
 
  L = assemLoad_LFE(Mesh,P1O2(),F,0);
  per = 0;
  progress_bar(per)
  for i = 1:NSTEPS
    if(per < floor(100*i/NSTEPS))
      per = floor(100*i/NSTEPS);
      progress_bar(per);
    end
    
    % Incorporate Dirichlet boundary data  
      
    [u_new,FreeDofs] = assemDir_LFE(Mesh,-1,G,i*dt);
    
    % Solve the linear system
    
    rhs = dt*L+S*u_old;
    rhs = rhs - M*u_new;
    u_new(FreeDofs) = M(FreeDofs,FreeDofs)\rhs(FreeDofs);
    
    L2Norm(i+1) = L2Err_LFE(Mesh,u_new,P1O2(),DUMMY);
    H1SNorm(i+1) = H1SErr_LFE(Mesh,u_new,P1O2(),DUMMY);
    LInfNorm(i+1) = LInfErr_LFE(Mesh,u_new,DUMMY);
    
    % Update vectors
    
    u_old = u_new;
    
  end
  
  % Generate figures
  
  plot_LFE(u_new,Mesh);
  colorbar;
  title('{\bf }');
  print('-depsc','Euler_Blowup_1.eps');
  
  fig = figure('Name','Blow up');
  plot(0:dt:T,L2Norm,'r-', ...
       0:dt:T,H1SNorm,'b-', ...
       0:dt:T,LInfNorm,'g-');
  title('{\bf Norm of solution}');
  xlabel('{\bf Time}');
  ylabel('{\bf Norm [log]}');
  legend('L2 norm','H1 semi-norm','LInf norm', ...
         'Location','NorthWest');
  set(gca,'YScale','log');
  print('-depsc','Euler_Blowup_2.eps');
  
  % Clear Memory
  
  clear all;
  