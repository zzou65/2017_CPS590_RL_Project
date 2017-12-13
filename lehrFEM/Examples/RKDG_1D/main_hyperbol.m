% Run script for linear advection problem.

% Copyright 2006-2007 Patrick Meury Holger Heumann
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland
  
  % Initialize constants
  X0 =0;                 % Left end point of interval
  X1 =2;                 % Right end point of interval
  NPOINTS =10.*[3];   % Initial set of grid points
  P=[2];                 % ploynomial degree 
  CFL=0.9*1./(2*P+1);
% CFL=0.1;
  T =2;                  % Final time
  M = 0;
%   
   U0 = @(x,t)(mod(x,2)>0.5).*(x<1.5);                 % Initial data
   U_ex =@(x,t)(mod(0.5+t,X1)<mod(1.5+t,X1))*...
                      ((x>mod(0.5+t,X1)).*(x<mod(1.5+t,X1)))+...
               (mod(0.5+t,X1)>mod(1.5+t,X1))*...
                      ((x>mod(0.5+t,X1)).*(x<X1)+(x>0).*(x<mod(1.5+t,X1)));      % Initial data
%     U0 = @(x,t)1+0.2*sin(pi*x);                % Initial data
%     U_ex =@(x,t)1+0.2*sin(pi*(x-t));           % Initial data
  G = @(x,t)zeros(size(x,1),1);                 % Right hand side load data
  FLUX = @(x) x;                                % Flux function
  NUMFLUX = @(x,y) x;                           % Numerical flux function
  ULIM = @mySlopeLim_lr;                          % Limiter for ansatz functions
  %ULIM = @(c,x,varargin)x;                      % no limiter
   
  % container for error calculation
  err=zeros(size(NPOINTS,2),size(P,2));
  h=err;
  
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
          %plot_hpDG_1D(Coordinates,p,M\L,@shap_Leg_1D);
  
          % Integrate ODE system
    
%         per = 0;
%         progress_bar(per)
  
          dx=min(Coordinates(2:end)-Coordinates(1:end-1));
          dt=CFL(q)*dx;
          NSTEPS=ceil(T/dt);
          Tend=NSTEPS*dt;
%           
          tau = zeros(1,NSTAGES+1);
          L = zeros(size(U,1),NSTAGES+1);
          V = zeros(size(U,1),NSTAGES+1);
          PV = zeros(size(U,1),NSTAGES+1);
          
          Tend2=0;
          for i = 1:NSTEPS 
             Tend2=Tend2+dt; 
%              if(per < floor(100*i/NSTEPS))
%                  per = floor(100*i/NSTEPS);
%                  progress_bar(per);
%              end      
%       
             % tau(1) = (i-1)*dt; 
             V(:,1) = U;
             
             PV(:,1) = ULIM(Coordinates, V(:,1),p, Inn_shap, Shap, QuadRule,M);
             %PV(:,1)=V(:,1);  
             
             Bvol=assemRKDG_update_Vol(Coordinates,PV(:,1),p,Shap,Grad_Shap,FLUX,QuadRule);
             Binn=assemRKDG_update_Inn(Coordinates,PV(:,1),p,Inn_shap,NUMFLUX);
             L(:,1) = M\(Bvol-Binn);
             
             for j = 2:(NSTAGES+1)
             %   tau(j) = 0;
                 V(:,j) = 0;

                 for k = 1:(j-1)
             %       tau(j) = tau(j) + A_RK(j-1,k)*tau(k) + dt*B_RK(j-1,k);
                     V(:,j) = V(:,j) + A_RK(j-1,k)*PV(:,k) + dt*B_RK(j-1,k)*(L(:,k));
                 end
      
                 if(j < NSTAGES+1)
                     %Lvol=assemLoad_hpDG_1D(Coordinates,p,QuadRule,Shap,G,tau(j));
                     
                     PV(:,j) = ULIM(Coordinates, V(:,j),p, Inn_shap, Shap, QuadRule,M);
                     %PV(:,j)=V(:,j);
                     
                     Bvol=assemRKDG_update_Vol(Coordinates,PV(:,j),p,Shap,Grad_Shap,FLUX,QuadRule);
                     Binn=assemRKDG_update_Inn(Coordinates,PV(:,j),p,Inn_shap,NUMFLUX);
                     L(:,j) = M\(Bvol-Binn);
                 end 
             end
                           
             % Update solution
             
             PV(:,NSTAGES+1) = ULIM(Coordinates, V(:,j),p, Inn_shap, Shap, QuadRule,M);
             %PV(:,NSTAGES+1)=V(:,NSTAGES+1);

             U = PV(:,NSTAGES+1);
%               TVDNorm_hpDG_1D(Coordinates,p,U,qr,err_Shap)
     
             if (mod(i,10)==0)
                 plot_hpDG_1D(Coordinates,p,U,@shap_Leg_1D);
                 set(gca, ... 
                     'XLim',[0,2], ...
                     'YLim',[-0.5,1.5]);
%                  L = assemLoad_hpDG_1D(Coordinates,p,QuadRule,Shap,U_ex,i*dt);
%                  plot_hpDG_1D(Coordinates,p,M\L,@shap_Leg_1D);
%                  set(gca, ... 
%                      'XLim',[0,2], ...
%                      'YLim',[-0.5,1.5]);
  
            end
            
          end          
          plot_hpDG_1D(Coordinates,p,U,@shap_Leg_1D);
          err(s,q)=L2Err_hpDG_1D(Coordinates,p,U,qr,err_Shap,U_ex,Tend);
        %  TVDNorm_hpDG_1D(Coordinates,p,U,qr,err_Shap)
          h(s,q)=dx;
      end % end of meshsize iteration
  end % end of polynomial degree iteration
  
  
  fig = figure('Name','Discretization error');
  hold on;
  for i=1:size(P,2)
    plot(h(:,i),err(:,i),'r+-');
  end
  hold off;
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  title('{\bf Discretization errors with respect to L^2 norm}');
  xlabel('{\bf Mesh width}');
  ylabel('{\bf Discretization error}');
  
   add_Slope(gca,'SouthEast',1);
   add_Slope(gca,'SouthEast',2);
   add_Slope(gca,'SouthEast',3);
  clear all;