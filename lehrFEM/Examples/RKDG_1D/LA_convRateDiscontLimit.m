% Run script for linear advection problem convergence plot
% Copyright 2006-2007 Patrick Meury Holger Heumann
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland
    
  % Initialize constants
  X0 =-0.5;                 % Left end point of interval
  X1 =2;                 % Right end point of interval
  NPOINTS =4*2.*[1:4];   % Initial set of grid points
  P=[0:2];                 % ploynomial degree 
  CFL=0.9*1./(2*P+1);
  T =0.5;                  % Final time  
  Mmm=10;                   % parameter for modified minmod   
%    U0 = @(x,t)(mod(x,2)>0.5).*(x<1.5);                 % Initial data
%    U_ex =@(x,t)(mod(0.5+t,X1)<mod(1.5+t,X1))*...
%                       ((x>mod(0.5+t,X1)).*(x<mod(1.5+t,X1)))+...
%                (mod(0.5+t,X1)>mod(1.5+t,X1))*...
%                       ((x>mod(0.5+t,X1)).*(x<X1)+(x>0).*(x<mod(1.5+t,X1)));      % Initial data
     U0 = @(x,t)((x>=0).*(x<=1));                                % Initial data
     U_ex=@(x,t)(((x-t)>=0).*((x-t)<=1));                        % Initial data
  G = @(x,t)zeros(size(x,1),1);                 % Right hand side load data
  FLUX = @(x) x;                                % Flux function
  NUMFLUX = @(x,y) x;                           % Numerical flux function
  ULIM = @mySlopeLim_lr;                          % Limiter for ansatz functions
  %ULIM = @(c,x,varargin)x;                      % no limiter
  %SULIM=@(c,x,varargin)x;
  
  % container for error calculation
  err=zeros(size(NPOINTS,2),size(P,2));
  h=err;
  err1=err;
  err2=err;
  
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
      qr=gauleg(-1,1,2*(ceil(P(q)+3/2)));
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
  
          % Integrate ODE system
     
          dx=min(Coordinates(2:end)-Coordinates(1:end-1));
          dt=CFL(q)*dx;
          NSTEPS=ceil(T/dt);
          Tend=NSTEPS*dt;
           
          tau = zeros(1,NSTAGES+1);
          L = zeros(size(U,1),NSTAGES+1);
          V = zeros(size(U,1),NSTAGES+1);
          PV = zeros(size(U,1),NSTAGES+1);
          
          for i = 1:NSTEPS 
             
             V(:,1) = U;
             
             PV(:,1) = ULIM(Coordinates, V(:,1),p, Inn_shap, Shap, QuadRule,Mmm);
             
             Bvol=assemRKDG_update_Vol(Coordinates,PV(:,1),p,Shap,Grad_Shap,FLUX,QuadRule);
             Binn=assemRKDG_update_Inn(Coordinates,PV(:,1),p,Inn_shap,NUMFLUX);
             L(:,1) = M\(Bvol-Binn);
             
             for j = 2:(NSTAGES+1)
                 V(:,j) = 0;

                 for k = 1:(j-1)
                     V(:,j) = V(:,j) + A_RK(j-1,k)*PV(:,k) + dt*B_RK(j-1,k)*(L(:,k));
                 end
      
                 if(j < NSTAGES+1)
                     PV(:,j) = ULIM(Coordinates, V(:,j),p, Inn_shap, Shap, QuadRule,Mmm);
                     
                     Bvol=assemRKDG_update_Vol(Coordinates,PV(:,j),p,Shap,Grad_Shap,FLUX,QuadRule);
                     Binn=assemRKDG_update_Inn(Coordinates,PV(:,j),p,Inn_shap,NUMFLUX);
                     L(:,j) = M\(Bvol-Binn);
                 end 
             end
                           
             % Update solution
                     
             PV(:,NSTAGES+1) = ULIM(Coordinates, V(:,j),p, Inn_shap, Shap, QuadRule,Mmm);
          
             U = PV(:,NSTAGES+1);
            
          end
          %plot_hpDG_1D(Coordinates,p,U,@shap_Leg_1D);
          err(s,q)=L2Err_hpDG_1D(Coordinates,p,U,qr,err_Shap,U_ex,Tend)
          err1(s,q)=L1Err_hpDG_1D(Coordinates,p,U,qr,err_Shap,U_ex,Tend);
          err2(s,q)=LinfErr_hpDG_1D(Coordinates,p,U,qr,err_Shap,U_ex,Tend);
          h(s,q)=dx;
      end % end of meshsize iteration
  end % end of polynomial degree iteration
  
  s=size(NPOINTS,2);
  fig = figure('Name','Discretization L^2-error,');
  plot(NPOINTS,err(:,1),'r+-',...
       NPOINTS,err(:,2),'b+-',...
       NPOINTS,err(:,3),'g+-')
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  xlabel('{\bf Mesh width}');
  ylabel('{\bf Discretization error}');
  legend('PolyDegree=0','PolyDegree=1','PolyDegree=2','Location','SouthEast');
  p = polyfit(log(NPOINTS(s-3:s)'),log(err(s-3:s,1)),1);
  add_Slope(gca,'NorthEast',p(1));
  p = polyfit(log(NPOINTS(s-3:s)'),log(err(s-3:s,2)),1);
  add_Slope(gca,'East',p(1));
  p = polyfit(log(NPOINTS(s-3:s)'),log(err(s-3:s,3)),1);
  add_Slope(gca,'SouthEast',p(1));
  
  
  print('-depsc', 'rate_L2_LADiscontLimit.eps');
  
  fig = figure('Name','Discretization L^1-error,');
  plot(NPOINTS,err1(:,1),'r+-',...
       NPOINTS,err1(:,2),'b+-',...
       NPOINTS,err1(:,3),'g+-')
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  xlabel('{\bf Mesh width}');
  ylabel('{\bf Discretization error}');
  legend('PolyDegree=0','PolyDegree=1','PolyDegree=2','Location','SouthEast');
  p = polyfit(log(NPOINTS(s-3:s)'),log(err1(s-3:s,1)),1);
  add_Slope(gca,'NorthEast',p(1));
  p = polyfit(log(NPOINTS(s-3:s)'),log(err1(s-3:s,2)),1);
  add_Slope(gca,'East',p(1));
  p = polyfit(log(NPOINTS(s-3:s)'),log(err1(s-3:s,3)),1);
  add_Slope(gca,'SouthEast',p(1));
  
  
  print('-depsc', 'rate_L1_LADiscontLimit.eps');
  
  fig = figure('Name','Discretization L^{\infty}-error,');
  plot(NPOINTS,err2(:,1),'r+-',...
       NPOINTS,err2(:,2),'b+-',...
       NPOINTS,err2(:,3),'g+-')
  grid('on');
  set(gca,'XScale','log','YScale','log','XDir','reverse');
  xlabel('{\bf Mesh width}');
  ylabel('{\bf Discretization error}');
  legend('PolyDegree=0','PolyDegree=1','PolyDegree=2','Location','SouthEast');
  p = polyfit(log(NPOINTS(s-3:s)'),log(err2(s-3:s,1)),1);
  add_Slope(gca,'NorthEast',p(1));
  p = polyfit(log(NPOINTS(s-3:s)'),log(err2(s-3:s,2)),1);
  add_Slope(gca,'East',p(1));
  p = polyfit(log(NPOINTS(s-3:s)'),log(err2(s-3:s,3)),1);
  add_Slope(gca,'SouthEast',p(1));
  
  print('-depsc', 'rate_Linf_LADiscontLimit.eps');
  
  clear all;