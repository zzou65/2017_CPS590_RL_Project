% initialize constants
%clear all

% Data I
A0_Handle=@(x,t,varargin) [(1-x(:,2)).*(1+x(:,2)).*sin(pi*x(:,1)) ...
                               (1-x(:,1)).*(1+x(:,1)).*sin(pi*x(:,2))];
V_Handle=@(x,varargin)[ones(size(x,1),1) 0.5.*ones(size(x,1),1)];                          

% Data II                               
U0_HANDLE=@(x,varargin)[sin(pi.*x(:,2)).*sin(pi.*x(:,1)) ...
    sin(pi.*x(:,1)).*sin(pi.*x(:,2))];
G_HANDLE=@(x,t,varargin)cos(2*pi*t).*[sin(pi.*x(:,2)).*sin(pi.*x(:,1)) ...
    sin(pi.*x(:,1)).*sin(pi.*x(:,2))];
F_HANDLE=@(x,t,varargin)cos(2*pi*t).*(...
    [+pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1))-pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2))...
     -pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1))+pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2))]+...
    [pi.^2.*cos(pi.*x(:,2)).*cos(pi.*x(:,1))+pi.^2.*sin(pi.*x(:,1)).*sin(pi.*x(:,2)) ...
     pi.^2.*cos(pi.*x(:,2)).*cos(pi.*x(:,1))+pi.^2.*sin(pi.*x(:,1)).*sin(pi.*x(:,2))])+...
    -2*pi*sin(2*pi*t).*(...
    [sin(pi.*x(:,2)).*sin(pi.*x(:,1)) ...
    sin(pi.*x(:,1)).*sin(pi.*x(:,2))]);
%     
VxCurl_Handle=@(x,t,varargin)cos(2*pi*t).*...
    [pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1))-pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2)) ...
    -pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1))+pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2))];    % - v x curl u

%Data III
% G_HANDLE = @(x,t,varargin)cos(2*pi*t)*[1-x(:,1) 1-x(:,2)];        % Dirichlet boundary data
% F_HANDLE = @(x,t,varargin)-2*pi*sin(2*pi*t)*[1-x(:,1) 1-x(:,2)];  % Right hand side source term        
% U0_HANDLE = @(x,varargin)[1-x(:,1) 1-x(:,2)];                     % Initial data

NREFS =4;
NSTEPS = 20;      % Number of time steps
NFRAMES = 200;    % Number of frames
NROUNDS = 2;      % Number of rounds
T = 1 ;            % Final time
THETA = 0.5;      % Coefficient of theta scheme
QuadRule=Duffy(TProd(gauleg(0,1,10)));
MU_HANDLE=@(x,varargin)1;
SIGMA_HANDLE = @(x,varargin)ones(size(x,1),1);
V_Handle=@(x,varargin)[ones(size(x,1),1) ones(size(x,1),1)];

% Open up FIFO buffer
 
buf = open(buffer(),'/scratch/users/hheumann');

Mesh.Coordinates = [-1 -1;1 -1;1 1;-1 1];
Mesh.Elements = [1 2 4;2 3 4];
Mesh = add_Edges(Mesh);
Mesh = add_Edge2Elem(Mesh);
    
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
Mesh.BdFlags(Loc) = [-1 -1 -2 -2];  % -1 inflowboundary, -2 noninflowboundary 
Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);  

for i = 1:NREFS
       Mesh=refine_REG(Mesh);
       Mesh=add_Edge2Elem(Mesh);         
end
Mesh=add_VBdFlags(Mesh);
Mesh=add_Patches(Mesh);


% precompute matrices
sigma=SIGMA_HANDLE(Mesh.Coordinates);
v=V_Handle(Mesh.Coordinates);

% sigma-weighted Mass of one forms
M_sigma=assemMat_W1Fdis(Mesh,sigma);

% topological curl operator
TopRot=assemMat_TopRot(Mesh);

% contraction of two forms 
ContrTwo=assemMat_ContrTwo(Mesh,v);

% diagonal mass matrix of two-forms
MassTwo=assemMat_MassTwoD(Mesh);

dt = T/NSTEPS;
S1=M_sigma+...
    dt*THETA*TopRot'*MassTwo*TopRot+...
    dt*M_sigma*ContrTwo*TopRot;
S2=M_sigma-...
    dt*(1-THETA)*TopRot'*MassTwo*TopRot-...
    dt*(1-THETA)*M_sigma*ContrTwo*TopRot;

 % Compute initial data
  
A_old=assemCochain_1f(Mesh,U0_HANDLE,gauleg(0,1,10));
     A_plot=plot_Norm_W1F_mov(A_old,Mesh);
     buf = push(buf,A_plot);
cL_old=assemCochain_1f(Mesh,F_HANDLE,gauleg(0,1,10),0);

% L2 error in every timestep;
err=zeros(1,NSTEPS+1);

per = 0;
progress_bar(per);
for i = 1:NSTEPS
    if(per < floor(100*i/NSTEPS))
      per = floor(100*i/NSTEPS);
      progress_bar(per);  
    end
    % L2 error
    err(i)=L2Err_W1F(Mesh,A_old,QuadRule,G_HANDLE,(i+1)*dt);
    
    % Assemble load vector 
    
    cL_new=assemCochain_1f(Mesh,F_HANDLE,gauleg(0,1,10),i*dt);
    
    % Incorporate boundary data  
    FreeDofs=intersect( find(Mesh.BdFlags ~=-1), find(Mesh.BdFlags ~=-2));
    nonInFlow=find(Mesh.BdFlags ~=-1);
    A_new=assemCochain_1f(Mesh,G_HANDLE,gauleg(0,1,10),i*dt);
    A_new(FreeDofs)=0;
    %!!! could be done faster
    cVxCurlBd=assemCochain_1f(Mesh,VxCurl_Handle,gauleg(0,1,10),i*dt);
    cVxCurlBd(nonInFlow)=0;
    
    % Solve the linear system
    
    rhs = S2*A_old + M_sigma*dt*THETA*cL_new+M_sigma*dt*(1-THETA)*cL_old;
    rhs = rhs - S1*A_new-dt*M_sigma*cVxCurlBd;
    A_new(FreeDofs) = S1(FreeDofs,FreeDofs)\rhs(FreeDofs);
    A_plot=plot_Norm_W1F_mov(A_new,Mesh);
    buf = push(buf,A_plot);

    %  plot_Norm_W1F(A_new,Mesh);
    %  colorbar;
    % %plot exact solution
    %  A_ex=assemCochain_1f(Mesh,G_Handle,gauleg(0,1,10),t*dt);
    %  A_ex_plot=plot_Norm_W1F_mov(A_ex,Mesh);
    %  buf2 = push(buf2,A_ex_plot);
    
    % coupling parameter
    
    % elctric field
    E=(-A_new+A_old)/dt;

    % electric current
    J=TopRot'*MassTwo*TopRot*A_new;

    % magnetic field
    B=TopRot*A_new;

    % Energy
    EJ=assemWedge_1f1f(Mesh,E,J);

    % ?? Momentum??
    BJ=assemWedge_2f1f(Mesh,B,J);
    
    % Update vectors
    
    A_old = A_new;
    cL_old = cL_new;
    
end

% L2 error
err(NSTEPS+1)=L2Err_W1F(Mesh,A_old,QuadRule,G_HANDLE,(i+1)*dt);

  fig = figure('Name','Error in timesteps');
  plot(0:dt:T,err,'ro--'); grid('on');
  set(gca,'YScale','log');
  xlabel('{\bf t}');
  ylabel('{\bf Error}');

  f=movie_LFE(buf,Mesh,NFRAMES,NROUNDS);
  close(buf);
  save movie f;
 clear all
