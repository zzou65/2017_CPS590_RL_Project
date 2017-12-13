% Run script for circular advection problem.

% Copyright 2006-2006 Patrick Meury
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland
  
  % Initialize constants
   
  NREFS = 4;             % Number of red refinement steps
  NSTEPS = 2000;         % Number of time steps
  AVIFILE = 'Blob.avi';  % Filename of .avi file           
  NFRAMES = 200;         % Number of movie frames
  T = 1;                 % Final time
  
  U0 = @u0_CA;           % Initial data
  UD = @uD_CA;           % Dirichlet boundary data
  G = @g_CA;             % Right hand side load data                           
  FLUX = @f_CA;          % Flux function
  
  NUMFLUX = @Upwind_CA;  % Numerical flux function
  ULIM = @limit_ad;      % Limiter for ansatz functions
  ALPHA = 4;             % Limiting parameter for adaptive limiter
  VLIM = 1;              % Limiter for test functions
  
  % Initialize TVD scheme (2-stage Runge-Kutta scheme)
    
  NSTAGES = 2;
  A_RK = [1 0; 1/2 1/2];
  B_RK = [1 0; 0   1/2];
  
  % Open FIFO buffer
  
  buf = open(buffer());
  
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
  CMin = min(U);
  CMax = max(U);
  buf = push(buf,U);
   
  % Integrate ODE system
    
  per = 0;
  progress_bar(per)    
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
      error('Time step scheme unstable - increase number of time steps');  
    end
    
    if(rem(i,max(floor(NSTEPS/NFRAMES),1)) == 0)
      buf = push(buf,U);
      CMax = max(CMax,max(U));
      CMin = min(CMin,min(U));
    end
    
  end
  clear('M','V','PV','L','Lvol','Bvol','Binn','Bbnd');
  
  % Generate auxiliary mesh
  
  nCoordinates = size(Mesh.Coordinates,1);
  nElements = size(Mesh.Elements,1);
  Coordinates = zeros(nCoordinates,2);
  Elements = zeros(nElements,3);
  V = zeros(nCoordinates,1);
  
  for i = 1:nElements
    Elements(i,:) = 3*(i-1)+[1 2 3];
    Coordinates(Elements(i,1),:) = Mesh.Coordinates(Mesh.Elements(i,1),:);
    Coordinates(Elements(i,2),:) = Mesh.Coordinates(Mesh.Elements(i,2),:);
    Coordinates(Elements(i,3),:) = Mesh.Coordinates(Mesh.Elements(i,3),:);
  end
  
  % Generate movie
  
  if(~isempty(AVIFILE))
    Mov = avifile(AVIFILE);
  end
  
  fig = figure('Name','Crouzeix-Raviart finite elements');
  XMax = max(Mesh.Coordinates(:,1));
  XMin = min(Mesh.Coordinates(:,1));
  XLim = [XMin XMax] + 0.05*(XMax-XMin)*[-1 1];
  YMax = max(Mesh.Coordinates(:,2));
  YMin = min(Mesh.Coordinates(:,2));
  YLim = [YMin YMax] + 0.05*(YMax-YMin)*[-1 1];
  CLim = [CMin CMax] + 0.05*(CMax-CMin)*[-1 1];
  
  pos = 1;
  while(~isempty(buf))
     
    % Extract solution from FIFO buffer  
       
    [U,buf] = pop(buf);
     
    V(Elements(:,1)) = U(Elements(:,2))+U(Elements(:,3))-U(Elements(:,1));
    V(Elements(:,2)) = U(Elements(:,1))+U(Elements(:,3))-U(Elements(:,2));
    V(Elements(:,3)) = U(Elements(:,1))+U(Elements(:,2))-U(Elements(:,3));
    
    % Generate figure
      
    clf;
    patch('Faces', Elements, ...
          'Vertices', [Coordinates(:,1) Coordinates(:,2) V], ...
          'CData', V, ...
          'FaceColor', 'interp', ...
          'EdgeColor', 'none');
    set(gca,'XLim',XLim, ...
            'YLim',YLim, ...
            'CLim',CLim, ...
            'DataAspectRatio',[1 1 1]);
    colorbar;
    drawnow();
    
    % Add frame to movie
    
    if(~isempty(AVIFILE))
      try
        F = getframe();
        Mov = addframe(Mov,F);
      catch
        fprintf('Invalid frame at position %d\n',pos);    
      end
    end
    
    pos = pos+1;
    
   end
   buf = close(buf);
    
   % Close the movie
  
  if(~isempty(AVIFILE))
    Mov = close(Mov);
  end
    
  % Clear memory
  
  clear all;
  