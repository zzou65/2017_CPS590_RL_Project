% 1D wave equation with absorbing boundary condition (ABC).
%
% Copyright 2006-2006 Patrick Meury
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize constants
  
  L = 2;                                                      % Length of the interval
  T = 1.5;                                                    % Final time
  C = 2;                                                      % Speed of sound
  ALPHA = 0.75;                                               % Right travelling wave
  BETA = 0.25;                                                % Left travelling wave
  U0 = @(x,varargin)(ALPHA+BETA)*pulse_1D(x,varargin{:});     % Initial data
  V0 = @(x,varargin)(ALPHA-BETA)*C*dpulse_1D(x,varargin{:});  % Initial velocity
  NPTS = 1000;                                                % Number of points
  NSTEPS = 10000;                                             % Number of time steps
   
  NFRAMES = 200;                                              % Number of frames for the movie
  FILENAME = 'Wave_1D_ABC.avi';                               % Filename of the .avi file
  XLim = [-L L];                                              % X-axes limits
  YLim = [-0.05 1.05];                                        % Y-axes limits
     
  % Initialize mesh
  
  dx = (XLim(2)-XLim(1))/(NPTS-1);
  Coordinates = transpose(XLim(1):dx:XLim(2));
  
  % Precompute matrices
  
  dt = T/NSTEPS;
  
  M = assemMat_P1_1D(Coordinates,@MASS_P1_1D);
  A = assemMat_P1_1D(Coordinates,@STIMA_Lapl_P1_1D);
  B = sparse([1 NPTS],[1 NPTS],[1 1]);

  S3 = 1/(dt)^2*M + C^2/4*A + C/(2*dt)*B;
  S2 = 2/(dt)^2*M - C^2/2*A;
  S1 = C/(2*dt)*B - 1/(dt)^2*M - C^2/4*A;
  
  % Open up FIFO buffer
 
  buf = open(buffer(),'/home/hheumann');
  
  % Compute initial data
  
  QuadRule = gauleg(0,1,2);
  
  rhs = assemLoad_P1_1D(Coordinates,QuadRule,U0);
  U1 = M\rhs;
  buf = push(buf,U1);
  
  % Compute initial velocity
  
  rhs = dt*assemLoad_P1_1D(Coordinates,QuadRule,V0)+M*U1;
  U2 = M\rhs;
  buf = push(buf,U2);
  
  % Integrate ODE system (Crank-Nicolson scheme, unconditionally stable)
  
  for i = 1:NSTEPS
      
    % Compute right hand side of the linear system
    
    rhs = S2*U2 + S1*U1;
    
    % Solve the linear system
    
    U3 = S3\rhs;
    if(rem(i,max(floor(NSTEPS/NFRAMES),1)) == 0)
      buf = push(buf,U3);
    end
    
    % Update solution vectors
    
    U1 = U2;
    U2 = U3;
      
  end
  
  % Generate movie
  
  if(~isempty(FILENAME))
    Mov = avifile(FILENAME);
  end
  fig = figure('Name','1D Wave equation with ABC');
  while(~isempty(buf))
    
    % Extract data from FIFO buffer  
      
    [U,buf] = pop(buf);
    
    % Generate figure
    
    clf;
    plot(Coordinates,U,'r-');
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
    
  % Clear memory and close figure
  
  clear all;
  close all;
  