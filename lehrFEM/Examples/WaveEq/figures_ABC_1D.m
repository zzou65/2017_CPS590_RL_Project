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
  
  FRAMES = [1 1500 3000 4500 6000 7500];                      % Frames of the plot
  FILENAME_1 = 'Wave_1D_ABC.eps';                             % Filename of the .eps file (wave plot)
  FILENAME_2 = 'Energy_1D_ABC.eps';                           % Filename of the .eps file (energy plot)      
  XLim = L*[-1 1];                                            % X-axes limits
  YLim = [-0.05 1.05*size(FRAMES,2)];                         % Y-axes limits
  
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
  
  % Preallocate memory
  
  nFrames = size(FRAMES,2);
  t = zeros(1,nFrames);
  U = zeros(NPTS,nFrames);
  Epot = zeros(1,NSTEPS);
  Ekin = zeros(1,NSTEPS);
  
  % Compute initial data
  
  QuadRule = gauleg(0,1,2);
  
  rhs = assemLoad_P1_1D(Coordinates,QuadRule,U0);
  U1 = M\rhs;
  frame = 1;
  if(~isempty(find(1 == FRAMES)))
    U(:,frame) = U1;
    t(frame) = 0;
    frame = frame+1;
  end
  
  % Compute initial velocity
  
  rhs = dt*assemLoad_P1_1D(Coordinates,QuadRule,V0)+M*U1;
  U2 = M\rhs;
  if(~isempty(find(2 == FRAMES)))
    U(:,frame) = U2;
    t(frame) = dt;
    frame = frame+1;
  end
  
  W = (U1+U2)/2;
  Epot(1) = C^2/2*transpose(W)*A*W;
  W = (U2-U1)/dt;
  Ekin(1) = 1/2*transpose(W)*M*W;
  
  % Integrate ODE system (Crank-Nicolson scheme, unconditionally stable)
  
  for i = 3:NSTEPS
      
    % Compute right hand side of the linear system
    
    rhs = S2*U2 + S1*U1;
    
    % Solve the linear system
    
    U3 = S3\rhs;
    if(~isempty(find(i == FRAMES)))
      U(:,frame) = U3;
      t(frame) = (i-1)*dt;
      frame = frame+1;
    end
    
    % Update solution vectors
    
    U1 = U2;
    U2 = U3;
     
    W = (U1+U2)/2;
    Epot(i-1) = C^2/2*transpose(W)*A*W;
    W = (U2-U1)/dt;
    Ekin(i-1) = 1/2*transpose(W)*M*W;
    
  end
  
  % Generate plot
  
  fig = figure('Name','1D Wave equation with ABC');
  hold on;
  for i = 1:nFrames
    plot(Coordinates,U(:,i)+(i-1),'r-');  
  end
  hold off;
  title('{\bf 1D Wave equation with ABC}');
  xlabel('{\bf x}');
  ylabel('{\bf t}');
  set(gca, ...
      'XLim',XLim, ...
      'YLim',YLim, ...
      'YTick',[], ...
      'YTickLabel',[]);
  print('-depsc',FILENAME_1);
  
  t = (1:NSTEPS)*dt;
  fig = figure('Name','Energy for 1D Wave equation with ABC');
  plot(t,Epot,'r-', ...
       t,Ekin,'b-', ...
       t,Ekin+Epot,'g-');
  legend('E_{pot}','E_{kin}','E_{tot}', ...
         'Location','NorthEast');
  title('{\bf Potential and kinetic energy}');
  xlabel('{\bf t}');
  ylabel('{\bf Energy}');
  set(gca, ...
      'XLim',[-0.05 1.05*T]);
  print('-depsc',FILENAME_2);
  
  % Clear memory and close figure
  
  clear all;
  