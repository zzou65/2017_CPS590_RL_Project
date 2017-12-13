function a=theta_error(A,L,T,U0,NPTS,NSTEPS,theta)

% 1D Wave equation with theta scheme
%
% Copyright 2008-2008 Holger Heumann
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

% Input parameters
%  
%   A=0                                                      % interval beginning 
%   L = 1;                                                    % Length of the interval
%   T = 1;                                                    % Final time
%   U0 = @(x,varargin)x.*(1-x);                       % Initial data
%   NPTS = 100;                                           % Number of points
%   NSTEPS = 100;                                        % Number of time steps
   
%   theta=1/2;
  
  % Initialize mesh
 
  dx = (L)/(NPTS-1);
  Coordinates = transpose(A:dx:A+L);
  
  % Precompute matrices
  
  dt = T/(NSTEPS-1); 
  QuadRule = gauleg(0,1,4);
  
  M = assemMat_P1_1D(Coordinates,@MASS_P1_1D);
  A = -assemMat_P1_1D(Coordinates,@STIMA_Lapl_P1_1D);
 
  S1 = M-dt^2*theta^2*A;
  S2 = M+dt^2*(1-theta)*theta*A;
  
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
  
    % Compute new value for V
   
    V_new = V_old+dt*(theta*A*U_new+(1-theta)*A*U_old);

    % Update old values for U and V
    
    U_old = U_new;
    V_old = V_new;
    
  end
  
  a=L2Err_P1_1D(Coordinates,U_new,gauleg(0,1,40),@(x,T,g)wave_solution(x,T,g),T,U0);
  return;