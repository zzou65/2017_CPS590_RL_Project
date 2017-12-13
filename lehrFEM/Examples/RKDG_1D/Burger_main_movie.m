% Run script for burger problem convergence plot

% Copyright 2006-2007 Patrick Meury Holger Heumann
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland
    
  % Initialize constants
  X0 =-5;                 % Left end point of interval
  X1 =10;                 % Right end point of interval
  NPOINTS =10.*10;    % Initial set of grid points
  P=0;               % ploynomial degree 
  CFL=0.9*1./(2*P+1);
% CFL=0.1;
  T =7;                  % Final time
  
  NFRAMES = 200;                                              % Number of frames
  FILENAME = 'StairB_1D_PML.avi';                               % Filename of .avi file
  XLim = [-5 10];                                    % X-axes limits
  YLim = [-0.2 1.2];   
   
  U0 = @(x,t)(x>0).*(x<1);                 % Initial data
  U_ex =@(x,t)((0<=t)*(t<=2)).*...
               ((x/t.*(0<x).*(x<t))+((t<x).*(x<1+t/2)))+...
               (t>2).*(x/t.*(0<x).*(x<(2*t)^(1/2)))    % Initial data
%       U0 = @(x,t)1+0.2*sin(pi*x);                % Initial data
%       U_ex =@(x,t)1+0.2*sin(pi*(x-t));           % Initial data
  G = @(x,t)zeros(size(x,1),1);                 % Right hand side load data
  FLUX = @(x) x.^2/2;                                % Flux function
  NUMFLUX = @(x,y) x.^2/2;                           % Numerical flux function
  %ULIM = @SlopeLim_lr;                          % Limiter for ansatz functions
 ULIM = @(c,x,varargin)x;                      % no limiter
  %SULIM=@(c,x,varargin)x;
  
  % Open up FIFO buffer
 
  buf = open(buffer(),'/home/hheumann/');

  % quadrature rule and shape functions at quadrature point and boundary
  for q=1:size(P,2)
      
      % Initialize TVD scheme (P+1-stage Runge-Kutta scheme)
 
      NSTAGES = min([P(q)+1 3]);
      %NSTAGES=3;
      A_RK=RKDG_alpha(NSTAGES);
      B_RK=RKDG_beta(NSTAGES);
        
      % quadrature rule and shape functions at quadrature point and boundary
      QuadRule=gauleg(-1,1,ceil(P(q)+3/2));
      Shap=shap_Leg_1D(QuadRule.x,P(q));
      Grad_Shap=grad_shap_Leg_1D(QuadRule.x,P(q));
      Inn_shap=shap_Leg_1D([-1;1],P(q));
      Inn_gradshap=grad_shap_Leg_1D([-1;1],P(q));
  
      % error calculation
      qr=gauleg(-1,1,100);
      err_Shap=shap_Leg_1D(qr.x,P(q));
  
      for s=1:size(NPOINTS,2);
      
          % Initialize mesh
      
          Coordinates = X0 + (X1-X0)/(NPOINTS(s)-1)*(0:(NPOINTS(s)-1));

          % Initialize polynomial degrees

          p=P(q)*ones(1,NPOINTS(s)-1);
  
          %Assemble mass matrix and extract diagonal
          M=assemMat_Vol_hpDG_1D(Coordinates,p,@MASS_Vol_hpDG_1D,QuadRule,Shap);
          %M=spdiags(M,0);
  
          % Compute initial data
  
          L = assemLoad_hpDG_1D(Coordinates,p,QuadRule,Shap,U0,0);
          U = M\L; 
          
          buf = push(buf,U);
  
          % Integrate ODE system
     
          dx=min(Coordinates(2:end)-Coordinates(1:end-1));
          dt=CFL(q)*dx;
          NSTEPS=ceil(T/dt);
          TVD=zeros(NSTEPS,1);
          Tend=NSTEPS*dt;
           
          tau = zeros(1,NSTAGES+1);
          L = zeros(size(U,1),NSTAGES+1);
          V = zeros(size(U,1),NSTAGES+1);
          PV = zeros(size(U,1),NSTAGES+1);
          for i = 1:NSTEPS 
             
             V(:,1) = U;
             
             PV(:,1) = ULIM(Coordinates, V(:,1),p, Inn_shap, Shap, QuadRule);
             
             Bvol=assemRKDG_update_Vol(Coordinates,PV(:,1),p,Shap,Grad_Shap,FLUX,QuadRule);
             Binn=assemRKDG_update_Inn(Coordinates,PV(:,1),p,Inn_shap,NUMFLUX);
             L(:,1) = M\(Bvol-Binn);
             
             for j = 2:(NSTAGES+1)
                 V(:,j) = 0;

                 for k = 1:(j-1)
                     V(:,j) = V(:,j) + A_RK(j-1,k)*PV(:,k) + dt*B_RK(j-1,k)*(L(:,k));
                 end
      
                 if(j < NSTAGES+1)
                     PV(:,j) = ULIM(Coordinates, V(:,j),p, Inn_shap, Shap, QuadRule);
                     
                     Bvol=assemRKDG_update_Vol(Coordinates,PV(:,j),p,Shap,Grad_Shap,FLUX,QuadRule);
                     Binn=assemRKDG_update_Inn(Coordinates,PV(:,j),p,Inn_shap,NUMFLUX);
                     L(:,j) = M\(Bvol-Binn);
                 end 
             end
                           
             % Update solution
                     
             PV(:,NSTAGES+1) = ULIM(Coordinates, V(:,j),p, Inn_shap, Shap, QuadRule);
          
             U = PV(:,NSTAGES+1);
             if(rem(i,max(floor(NSTEPS/NFRAMES),1)) == 0)
                buf = push(buf,U);
             end

             TVD(i)=TVDNorm_hpDG_1D(Coordinates,p,U,qr,err_Shap);
             if (mod(i,10)==0)
                 plot_hpDG_1D(Coordinates,p,U,@shap_Leg_1D);
                 set(gca, ... 
                     'XLim',[-5,10], ...
                     'YLim',[-0.2,1.2]); 
%                  L = assemLoad_hpDG_1D(Coordinates,p,QuadRule,Shap,U_ex,i*dt);
%                  plot_hpDG_1D(Coordinates,p,M\L,@shap_Leg_1D);
%                  set(gca, ... 
%                      'XLim',[-5,10], ...
%                      'YLim',[-0.5,1.5]);

             end
            
          end
          fig = figure('Name','TVDM-Norm,');
          plot(1:NSTEPS,TVD);
          title('{\bf Discretization errors with respect to L^2 norm}');
          xlabel('{\bf TimeStep}');
          ylabel('{\bf TVDM norm}');
  % Generate movie
  
  if(~isempty(FILENAME))
    Mov = avifile(FILENAME);
  end
  fig = figure('Name','advection');
  while(~isempty(buf))
    
     % Extract data from FIFO buffer  
      
     [U,buf] = pop(buf);
    
     % Generate figure
    
     clf;
      plot_hpDG_1D(Coordinates,p,U,@shap_Leg_1D);
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
  
  end % end of meshsize iteration
  end % end of polynomial degree iteration
  
  
  % Clear memory and closes figure
  
  clear all;