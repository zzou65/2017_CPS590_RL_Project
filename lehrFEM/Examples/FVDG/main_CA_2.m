% Run script for circular advection problem.

% Copyright 2006-2006 Patrick Meury
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland
  
  % Initialize constants
   
  NREFS = 5;                   % Number of red refinement steps
  NSTEPS = 2000;               % Number of time steps          
  FRAMES = [0 200:200:2000];   % Number of movie frames
  T = 1;                       % Final time
  
  U0 = @u0_CA;                 % Initial data
  UD = @uD_CA;                 % Dirichlet boundary data
  UEX = @uex_CA;               % Exact solution
  G = @g_CA;                   % Right hand side load data                           
  FLUX = @f_CA;                % Flux function
  
  NUMFLUX = @Upwind_CA;        % Numerical flux function
  %ULIM = @(Mesh,U,varargin)U;
  %ULIM = @limit_ad;            % Limiter for ansatz functions
  ULIM = @limit_p0;
  ALPHA = 4;                   % Limiting parameter for adaptive limiter
  %VLIM = 1;                    % Limiter for test functions
  VLIM = 0;
  
  % Initialize TVD scheme (2-stage Runge-Kutta scheme)
    
  NSTAGES = 2;
  A_RK = [1 0; 1/2 1/2];
  B_RK = [1 0; 0   1/2];
  
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_Circ.dat','Elem_Circ.dat');
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);         
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
  Mesh.BdFlags(Loc) = -1;
  for i = 1:NREFS
    Mesh = refine_REG(Mesh,@dist_circ,[0 0],1);  
  end
  Mesh = add_Edge2Elem(Mesh);
  Mesh = add_DGData(Mesh);
  
  % Assemble mass matrix
  
  QuadRule_2D = P3O3();
  QuadRule_1D = gauleg(0,1,2);
  
  M = diag(assemMat_Vol_DG(Mesh,@MASS_Vol_DGCR));
         
  % Compute initial data
  
  L = assemLoad_Vol_DG(Mesh,@LOAD_Vol_DGCR,QuadRule_2D,U0);
  U = L./M;
   
  if(~isempty(find(0 == FRAMES)))
      
    plot_DGCR(ULIM(Mesh,U,ALPHA),Mesh);
    title('{\bf Approximate solution}');
    xlabel(['{\bf ' int2str(0) ' time steps}']);
    set(gca,'XLim',[-1.1 1.1], ...
            'YLim',[-1.1 1.1], ...
            'ZLim',[-0.1 1.1], ...
            'CLim',[-0.1 1.1]);
    colorbar();
    drawnow();
    print('-depsc',['CAdv_' int2str(0) '.eps']);
    close all;
    drawnow();
    
    contour_DGCR(U,Mesh,.1:.1:.9,'colorbar');
    hold on;
    plot(cos(2*pi*(0:.01:1)),sin(2*pi*(0:.01:1)),'k-');
    hold off;
    title('{\bf Approximate solution}');
    xlabel(['{\bf ' int2str(0) ' time steps}']);
    set(gca,'XLim',[-1.1 1.1], ...
            'YLim',[-1.1 1.1]);
    drawnow();
    print('-depsc',['CAdv_LC_' int2str(0) '.eps']);
    close all;
    drawnow();
    
  end
  
  % Integrate ODE system
    
  per = 0;
  progress_bar(per);    
  dt = T/NSTEPS;
  tau = zeros(1,NSTAGES+1);
  L = zeros(size(U,1),NSTAGES+1);
  V = zeros(size(U,1),NSTAGES+1);
  PV = zeros(size(U,1),NSTAGES+1);
  L1Err = zeros(1,NSTEPS);
  L2Err = zeros(1,NSTEPS);
  for i = 1:NSTEPS 
    if(per < floor(100*i/NSTEPS))
      per = floor(100*i/NSTEPS);
      progress_bar(per);
    end  
       
    tau(1) = (i-1)*dt; 
    V(:,1) = U;
    
    Lvol = assemLoad_Vol_DG(Mesh,@LOAD_Vol_DGCR,QuadRule_2D,G,tau(1));
    
    PV(:,1) = ULIM(Mesh,V(:,1),ALPHA);
    Bvol = assemConv_Vol_DGCR(Mesh,PV(:,1),VLIM,QuadRule_2D,FLUX);
    Binn = assemConv_Inn_DGCR(Mesh,PV(:,1),VLIM,QuadRule_1D,NUMFLUX);
    Bbnd = assemConv_Bnd_DGCR(Mesh,PV(:,1),VLIM,QuadRule_1D,NUMFLUX,UD,tau(1));
    
    L(:,1) = Lvol+Binn+Bbnd+Bvol;
      
    for j = 2:(NSTAGES+1)
      tau(j) = 0;
      V(:,j) = 0;
 
      for k = 1:(j-1)
        tau(j) = tau(j) + A_RK(j-1,k)*tau(k) + dt*B_RK(j-1,k);
        V(:,j) = V(:,j) + A_RK(j-1,k)*PV(:,k) + dt*B_RK(j-1,k)*(L(:,k)./M);
      end
      
      if(j < NSTAGES+1)
        Lvol = assemLoad_Vol_DG(Mesh,@LOAD_Vol_DGCR,QuadRule_2D,G,tau(j));
        
        PV(:,j) = ULIM(Mesh,V(:,j),ALPHA);
        Bvol = assemConv_Vol_DGCR(Mesh,PV(:,j),VLIM,QuadRule_2D,FLUX);
        Binn = assemConv_Inn_DGCR(Mesh,PV(:,j),VLIM,QuadRule_1D,NUMFLUX);
        Bbnd = assemConv_Bnd_DGCR(Mesh,PV(:,j),VLIM,QuadRule_1D,NUMFLUX,UD,tau(j));
        
        L(:,j) = Lvol+Binn+Bbnd+Bvol;
      end
      
    end
         
    % Update solution
   
    U = V(:,NSTAGES+1);
    if(max(abs(U)) > 10)
      error('Time step scheme unstable - increase number of time steps');  
    end
    L1Err(i) = L1Err_DGCR(Mesh,U,QuadRule_2D,UEX,i*dt);
    L2Err(i) = L2Err_DGCR(Mesh,U,QuadRule_2D,UEX,i*dt);
    if(~isempty(find(i == FRAMES)))   
  
      plot_DGCR(U,Mesh);
      title('{\bf Approximate solution}');
      xlabel(['{\bf ' int2str(i) ' time steps}']);
      set(gca,'XLim',[-1.1 1.1], ...
              'YLim',[-1.1 1.1], ...
              'ZLim',[-0.1 1.1], ...
              'CLim',[-0.1 1.1]);
      colorbar();
      drawnow();
      print('-depsc',['CAdv_' int2str(i) '.eps']);
      close all;
      drawnow();
      
      contour_DGCR(U,Mesh,.1:.1:.9,'colorbar');
      hold on;
      plot(cos(2*pi*(0:.01:1)),sin(2*pi*(0:.01:1)),'k-');
      hold off;
      title('{\bf Approximate solution}');
      xlabel(['{\bf ' int2str(i) ' time steps}']);
      set(gca,'XLim',[-1.1 1.1], ...
              'YLim',[-1.1 1.1]);
      drawnow();
      print('-depsc',['CAdv_LC_' int2str(i) '.eps']);
      close all;
      drawnow();
      
    end
    
  end
  clear('M','V','PV','L','Lvol','Bvol','Binn','Bbnd');
       
  % Generate figure
  
  figure('Name','Discretization error');
  plot(1:NSTEPS,L1Err,'r-', ...
       1:NSTEPS,L2Err,'b-');
  legend('L1 error', ...
         'L2 error', ...
         'Location','SouthEast');
  title('{\bf Discretization errors}');
  xlabel('{\bf Time step}');
  ylabel('{\bf Discretization error [log]}');
  set(gca,'YScale','log');
  drawnow();
  print('-depsc','DiscErr.eps');
  close all;
  
  % Clear memory
  
  clear all;
  