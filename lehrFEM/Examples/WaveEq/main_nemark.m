% 1D Wave equation with perfectly matched layer (PML)
%
% Copyright 2006-2006 Patrick Meury
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize constants
  
  A=0                                                      % interval beginning 
  L = 1;                                                    % Length of the interval
  T = 1;                                                    % Final time
  C = 1;                                                      % Speed of sound
  U0 = @(x,varargin)x.*(1-x);     % Initial data
  NPTS = 100;                                                % Number of points
  NSTEPS = 100;                                             % Number of time steps
   
  NFRAMES = 200;                                              % Number of frames
  FILENAME = 'Wave_1D_PML.avi';                               % Filename of .avi file
  XLim = [A A+L];                                    % X-axes limits
  YLim = [-1.05 1.05];                                        % Y-axes limits
  
  gamma=1/2;
  beta=1;
  
  newmark_error(A,L,T,U0,NPTS,NSTEPS,beta,gamma)
  % Initialize mesh
 
  dx = (L)/(NPTS-1);
  Coordinates = transpose(A:dx:A+L);
  
  % Precompute matrices
  
  dt = T/(NSTEPS-1); 
  QuadRule = gauleg(0,1,4);
  
  M = assemMat_P1_1D(Coordinates,@MASS_P1_1D);
  A = -assemMat_P1_1D(Coordinates,@STIMA_Lapl_P1_1D);
 
  S1 = M-dt^2*beta*A;
  S2 = M+dt^2*(1/2-beta)*A;
  % Open up FIFO buffer
 
  buf = open(buffer(),'/home/hheumann/');

  % Compute initial data
 
  U_old = M\assemLoad_P1_1D(Coordinates,QuadRule,U0);
  V_old = zeros(NPTS,1);
  
  % Integrate ODE system (dissipative leapfrog scheme, CFL condition)
  
  FreeDofs = 2:(NPTS-1);
  for i = 2:NSTEPS
    
    % Compute right hand side load data
    
    rhs = S2*U_old + dt*V_old;
    
    % Incorporate zero Dirichlet boundary data
    
    U_new = zeros(NPTS,1);
    
    % Compute new value for U
    
    U_new(FreeDofs) = S1(FreeDofs,FreeDofs)\rhs(FreeDofs);
     if(rem(i,max(floor(NSTEPS/NFRAMES),1)) == 0)
      buf = push(buf,U_new);
    end
    % Compute new value for V
   
    V_new = V_old+dt*(gamma*A*U_new+(1-gamma)*A*U_old);

    % Update old values for U and V
    
    U_old = U_new;
    V_old = V_new;
    
  end
  
  % Generate movie
  
  if(~isempty(FILENAME))
    Mov = avifile(FILENAME);
  end
  fig = figure('Name','1D Wave equation with PML');
  while(~isempty(buf))
    
    % Extract data from FIFO buffer  
      
    [U,buf] = pop(buf);
    
    % Generate figure
    
    clf;
    plot(Coordinates,U,'r-', ...
         [-L -L],YLim,'k--', ...
         [ L  L],YLim,'k--');
    set(gca, ...
        'XLim',XLim, ...
        'YLim',YLim);
    
    % Add frame to movie
    
    if(~isempty(FILENAME))
      F = getframe();
      Mov = addframe(Mov,F);
    end
    
  end
    
  % Close the movie
  
  if(~isempty(FILENAME))
    Mov = close(Mov);
  end
  
  % Clear memory and closes figure
  
  clear all;
  close all;