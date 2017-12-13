% Run script for circular advection problem.

% Copyright 2006-2006 Patrick Meury
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland
  
  % Initialize constants
   
  NREFS = 5;             % Number of red refinement steps
  NSTEPS = 100;          % Number of time steps
  T = 0.25;              % Final time
  U0 = @u0_CA;           % Initial data
  UD = @uD_CA;           % Dirichlet boundary data
  UEX = @uex_CA;         % Exact solution
  G = @g_CA;             % Right hand side load data                           
  FLUX = @f_CA;          % Flux function
  NUMFLUX = @Upwind_CA;  % Numerical flux function
  TYPE = 1;              % Discretization type
  
  % Initialize limiters and time stepping schemes
  
  switch(TYPE)
    case 0    
      ULIM = @limit_p0;
      ALPHA = 0;
      VLIM = 0;
      NSTAGES = 1;
      A_RK = 1;
      B_RK = 1;
    case 1    
      ULIM = @limit_ad;
      ALPHA = 4;
      VLIM = 1;
      NSTAGES = 2;
      A_RK = [1 0; 1/2 1/2];
      B_RK = [1 0; 0   1/2];
  end
      
  % Initialize quadrature rules
  
  QuadRule_1D = gauleg(0,1,2);
  QuadRule_2D = P3O3();
  
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_Circ.dat','Elem_Circ.dat');
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  Mesh = add_Edges(Mesh);         
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
  Mesh.BdFlags(Loc) = -1;
  
  h = zeros(1,NREFS);
  L1Err = zeros(1,NREFS);
  L2Err = zeros(1,NREFS);
  for ref = 1:NREFS
  
    % Refine mesh and add DG data  
      
    Mesh = refine_REG(Mesh,@dist_circ,[0 0],1);  
    Mesh = add_Edge2Elem(Mesh);
    Mesh = add_DGData(Mesh);
  
    % Compute initial data
    
    M = diag(assemMat_Vol_DG(Mesh,@MASS_Vol_DGCR));
    L = assemLoad_Vol_DG(Mesh,@LOAD_Vol_DGCR,QuadRule_2D,U0);
    U = L./M;
  
    % Integrate ODE system
    
    per = 0;
    progress_bar(per);    
    dt = T/NSTEPS;
    tau = zeros(1,NSTAGES+1);
    L = zeros(size(U,1),NSTAGES+1);
    V = zeros(size(U,1),NSTAGES+1);
    PV = zeros(size(U,1),NSTAGES+1);
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
        error('Time stepping scheme unstable - increase number of time steps');  
      end
      
    end
    clear('M','V','PV','L','Lvol','Bvol','Binn','Bbnd');
  
    h(ref) = get_MeshWidth(Mesh);
    L1Err(ref) = L1Err_DGCR(Mesh,U,QuadRule_2D,UEX,T);
    L2Err(ref) = L2Err_DGCR(Mesh,U,QuadRule_2D,UEX,T);
  
    % Update time steps
    
    NSTEPS = 2*NSTEPS;
    
  end
  
  % Generate figure
  
  figure('Name','Discretization error');
  plot(h,L1Err,'r-', ...
       h,L2Err,'b-', ...
       h,L1Err,'k+', ...
       h,L2Err,'k+');
  legend('L1 error', ...
         'L2 error', ...
         'Location','SouthEast');
  title('{\bf Discretization errors}');
  xlabel('{\bf Mesh Width [log]}');
  ylabel('{\bf Discretization error [log]}');
  set(gca,'XScale','log','XDir','reverse','YScale','log');
  
  % Clear memory
  
  % clear all;
  