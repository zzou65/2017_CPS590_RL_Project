function [U,L2normU,LmaxnormU] = LagrangianMainSolver(mesh,Tfinal,Diff,N,Lhandle,Dhandle,Ihandle,Chandle,dosearch)
% LagrangianFinalSolver Computes a 2D Parabolic PDE of Convection Diffusion 
% form.
%
%   U = LagrangianFinalSolver(MESH,Tfinal,Diff,N,Lhandle,Dhandle,Ihandle,Chandle,dosearch) solves a 
%   2D convection Diffusion PDE for the given Mesh with Diffusion parameter Diff at the
%   Tfinal time instant with N divisions on Tfinal.
%
%
% Lhandle is for Source Function              e.g Lhandle = @(x,y,t,varargin) 0; 
% Dhandle is for Dirichlet Boundary Function  e.g Dhandle = @(x,varargin) 0;     
% Ihandle is for Initial Condition             
% e.g Ihandle = @(x,varargin) (max(0,1-4*sqrt((x(:,1)-0.5).^2+(x(:,2)-0.25).^2))); 
% Chandle is for Convection vector at a point
% e.g = @(x,varargin) [-sin(pi*x(1))*cos(pi*x(2)) sin(pi*x(2))*cos(pi*x(1))];

% dosearch = 0/1 if velocity is time independent then search is performed
% only once. If velocity is time independent dosearch = 0 else dosearch = 1

%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%
%
%   Example:
%
%   Dt(U) + Dx(U) + Dy(U) = Diff*Lapacian(U) + L(x,t);
%   U = assemLoad_LFE(Mesh,Tfinal,Diff,N);
%
%   By Rajdeep Deb
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
%Timestep Calculation
tau = (Tfinal/N);

%Dirichlet Boundary Points and Values Calculated
[UD,F] = assemDir_LFE(mesh,-1,Dhandle);

%Initial Condiations applied
U0 = zeros(size(mesh.Coordinates(:,1),1),1);
U0(:,1) = Dhandle(mesh.Coordinates(:,1),mesh.Coordinates(:,2));
U0(F) = Ihandle(mesh.Coordinates(F,:));

fig = figure('name','Semi-Lagrangian LFE','renderer','painters');
plot_LFE(U0,mesh,fig); colorbar; 
title('Initial Data'); drawnow;

%Lagrangian Mass matrix, Laplacian Matrix and Load Matrix Cooresponding to
%Previous time step and Lagrangian tracking position is calculated

[M,S,L,storeElement] = SemiLagU_LFE(mesh,U0,tau,Chandle,Dhandle);

U = U0;
UT = U0;
L2normU = zeros(1,N);
L2normU(1) = norm(U0,2);

LmaxnormU = zeros(1,N);
LmaxnormU(1) = max(abs(U0));
timestep = 1;

%Loop through each time step to obtaion the solution
for T = 0:tau:(Tfinal - tau)
  fprintf('timestep %d, t = %f\n',timestep,T+tau);
    timestep = timestep + 1;
    L2 = tau*SourceLoad(mesh,P7O6,Lhandle,T + tau);
    
    S = Diff*S;
    K = (M + S)\(M*L + L2);
    U(F) = K(F);
    
    UT(F) = U(F);
    if(dosearch == 1)
        [M,S,L] = SemiLagU_LFE(mesh,UT,tau,Chandle,Dhandle);
    else
        [L] = SemiLagUConstantVel_LFE(mesh,UT,tau,Chandle,Dhandle,storeElement);
    end

    plot_LFE(U,mesh,fig); colorbar; 
    title(sprintf('t = %f',T+tau)); drawnow; 
    
    L2normU(timestep) = norm(UT,2);
    LmaxnormU(timestep) = max(abs(UT));
end;