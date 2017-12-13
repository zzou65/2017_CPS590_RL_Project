% Test Case to solve a given problem using SemiLagrangian LFE 
mesh = Sqr_Mesh(4); % Design mesh on the domain
Tfinal = 1;         % Solution to be obtained at Tfinal
Diff = 0.0;           % Diffusion coefficient of the PDE
NumberTimeSteps = 20; % Number of time steps

Lhandle = @(x,y,t,varargin) 0; %Source Load
Dhandle = @(x,varargin) 0;     %Dirichlet Boundary Function
Ihandle = @(x,varargin) (max(0,1-4*sqrt((x(:,1)-0.5).^2+(x(:,2)-0.25).^2))); %Initial Condition Function
Chandle = @(x,varargin) [-sin(pi*x(1))*cos(pi*x(2)) sin(pi*x(2))*cos(pi*x(1))]; %Convection Functions

dosearch = 0; %Must be zero if velocity is time independent. Then only once search is performed

[U,L2normU,LmaxnormU,M,S,L] = LagrangianMainSolver(mesh,Tfinal,Diff,NumberTimeSteps,Lhandle,Dhandle,Ihandle,Chandle,dosearch);

figure(NumberTimeSteps+3);
plot(L2normU)
title('L2norm of solution with time');

figure(NumberTimeSteps+4);
plot(LmaxnormU)
title('Linfinite norm of solution with time');