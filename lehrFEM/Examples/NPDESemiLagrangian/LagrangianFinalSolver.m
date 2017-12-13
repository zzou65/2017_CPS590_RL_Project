function [U] = LagrangianFinalSolver(mesh,Tfinal,Diff,N)
% LagrangianFinalSolver Computes a 2D Parabolic PDE of Convection Diffusion 
% form.
%
%   U = LagrangianFinalSolver(MESH,Tfinal,Diff,N) solves a 2D convection 
%   diffusion PDE for the given Mesh with Diffusion parameter Diff at the
%   Tfinal time instant with N divisions on Tfinal.
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

%Initial, Boundary,Load and Velocity Vector Handler Creation
Lhandle = @(x,y,t,varargin) 0; %Source Load
Dhandle = @(x,varargin) 0;     %Dirichlet Boundary Function
Ihandle = @(x,varargin) exp(-(1/Diff)*((x(:,1) - 0.5).^2 + (x(:,2) - 0.5).^2)); %Initial Condition Function
Chandle = @(x,varargin) [1 1]; %Convection Functions

%Dirichlet Boundary Points and Values Calculated
[UD,F] = assemDir_LFE(mesh,-1,Dhandle);

%Initial Condiations applied
U0 = zeros(size(mesh.Coordinates(:,1),1),1);
U0(:,1) = Dhandle(mesh.Coordinates(:,1),mesh.Coordinates(:,2));
U0(F) = Ihandle(mesh.Coordinates(F,:));

plot_LFE(U0,mesh);colorbar

%Lagrangian Mass matrix, Laplacian Matrix and Load Matrix Cooresponding to
%Previous time step and Lagrangian tracking position is calculated
[M,S,L] = SemiLagrangian_LFE(mesh,U0,tau,Chandle,Dhandle);

U = U0;
UT = U0;

%Loop through each time step to obtaion the solution
for T = 0:tau:(Tfinal - tau)
    L2 = tau*SourceLoad(mesh,P7O6,Lhandle,T + tau);
    
    S = Diff*S;
    K = (M + S)\(L + L2);
    U(F) = K(F);
    
    UT(F) = U(F);

    [M,S,L] = SemiLagrangian_LFE(mesh,UT,tau,Chandle,Dhandle);
    plot_LFE(U,mesh);colorbar
end;

%Plot the result
plot_LFE(U,mesh);colorbar