function [] = example1()
%EXAMPLE1 plane wave DG example
%
%   Solves a homogeneous Helmholtz problem with impedance, Dirichlet and
%   Neumann boundary conditions.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland


  %% Generate mesh

  % Initialize mesh constants
  R = 1;
  BBOX = [-R -R; R R];                                                            % Bounding box
  H0 = 0.1;                                                                       % Initial mesh width
  DHANDLE = @(x) dist_diff(dist_circ(x,[0 0],R),...
    dist_union(dist_circ(x,[0 0.3],0.2),dist_rect(x,[-0.3 -0.4],0.6,0.2)));       % Signed distance function
%   DHANDLE = @(x) dist_diff(dist_rect(x,[-R -R],2*R,2*R),...
%     dist_union(dist_circ(x,[0 0.3],0.2),dist_rect(x,[-0.3 -0.4],0.6,0.2)));       % Signed distance function
  HHANDLE = @h_uniform;                                                           % Element size function
%   HHANDLE = @(x,varargin) max(1,sqrt(sum(x.^2,2)));
  FIXEDPOS = [-0.3 -0.4; 0.3 -0.4; 0.3 -0.2; -0.3 -0.2];                          % Fixed boundary vertices of the mesh
%   FIXEDPOS = [-0.3 -0.4; 0.3 -0.4; 0.3 -0.2; -0.3 -0.2; -R -R; R -R; R R; -R R];  % Fixed boundary vertices of the mesh
  DISP = 1;                                                                       % Display flag
  IMPBD = -1;
  DIRBD = -2;
  NEUBD = -3;
  
  % Construct mesh
  Mesh = init_Mesh(BBOX,H0,DHANDLE,HHANDLE,FIXEDPOS,DISP);  
%   Mesh = morph(Mesh);
  Mesh = add_Edges(Mesh);
  
  % Set boundary edges
  Loc = get_BdEdges(Mesh);
  pnt = Mesh.Coordinates(Mesh.Edges(Loc,1),:);
  indImp = sum(pnt.^2,2)>0.9;
  LocImp = Loc(indImp);
  Loc = Loc(~indImp);
  pnt = pnt(~indImp,:);
  LocDir = Loc(pnt(:,2)>=0);
  LocNeu = Loc(pnt(:,2)<0);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(LocImp) = IMPBD;
  Mesh.BdFlags(LocDir) = DIRBD;
  Mesh.BdFlags(LocNeu) = NEUBD;
  
  % Add other data
  Mesh = orient_Elems(Mesh);
  Mesh = add_Edge2Elem(Mesh);
  Mesh = add_DGData(Mesh);
  
  % Plot mesh
  plot_Mesh(Mesh,'as');
  
  %% Set problem parameters
  
  % Set boundary data
  omega = 3*pi;
  d = [2 -1];
  d = d/norm(d);
  gD = @(x,varargin) exp(i*omega*x*d');
  gI = @(x,n,varargin) zeros(size(x,1),1);
  gN = @(x,n,varargin) zeros(size(x,1),1);
  
  % Add plane wave data to mesh
  Mesh = set_Data_PWDG(Mesh,'nDofs',5,'MakeDir','rand0','Omega',omega,...
  'a',@(w,h,varargin) 2/(w*h),'b',@(w,h,varargin) 0.01*w*h);

  %% Assemble data
  
  % Define quadrature rule
  qr1 = gauleg(0,1,12);
  
  % Assemble stiffness matrix
  [I_Inn,J_Inn,A_Inn] = assemMat_Inn_PDG2(Mesh,@STIMA_Helm_Inn_PWDG,omega);
  [I_Imp,J_Imp,A_Imp] = assemMat_Bnd_PDG2(Mesh,IMPBD,@STIMA_Helm_Imp_Bnd_PWDG,omega);
  [I_Dir,J_Dir,A_Dir] = assemMat_Bnd_PDG2(Mesh,DIRBD,@STIMA_Helm_Dir_Bnd_PWDG,omega);
  [I_Neu,J_Neu,A_Neu] = assemMat_Bnd_PDG2(Mesh,NEUBD,@STIMA_Helm_Neu_Bnd_PWDG,omega);
  A = sparse([I_Inn;I_Imp;I_Dir;I_Neu],...
    [J_Inn;J_Imp;J_Dir;J_Neu],...
    [A_Inn;A_Imp;A_Dir;A_Neu]);
  
  % Assemble load vector
  b_Imp = assemLoad_Bnd_PDG2(Mesh,IMPBD,@LOAD_Imp_Bnd_PWDG,qr1,omega,gI);
  b_Dir = assemLoad_Bnd_PDG2(Mesh,DIRBD,@LOAD_Dir_Bnd_PWDG,qr1,omega,gD);
  b_Neu = assemLoad_Bnd_PDG2(Mesh,NEUBD,@LOAD_Neu_Bnd_PWDG,qr1,omega,gN);
  b = b_Imp + b_Dir + b_Neu;
  
  
  %% Solve system
  
  % Construct preconditioner
  [PL,PR] = assemPrec_SVD_PDG(Mesh,A);
  
  % Solve equation with preconditioning
  w = (PL*A*PR)\(PL*b);
  u = PR*w;
  
  % Plot solution
  plot_PWDG(u,Mesh,omega,2,[1 0 1]);
  plot_PWDG_Dir(Mesh,u);
  contour_PWDG(u,Mesh,[1 0 1]);
  