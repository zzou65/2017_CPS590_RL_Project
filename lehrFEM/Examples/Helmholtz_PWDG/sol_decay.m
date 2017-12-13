function [] = sol_decay()
%SOL_DECAY Check for decay in periodic solution
%   Detailed explanation goes here

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Define problem parameters
omega = 2*pi;
d = [-1 0];
d = d/norm(d);
psi = 0;
u_ex = @(x,varargin) exp(i*(omega*x*d'+psi));
% du_ex = @(x,varargin) i*omega*exp(i*(omega*x*d'+psi)).*d(ones(size(x,1),1),:);
g = @(x,n,varargin) i*omega*(1+d*n')*exp(i*(omega*x*d'+psi));
% f = @(x,varargin) zeros(size(x,1),1);

% Define discretization parameters
nref = 5;
ndir = 3;
phi = 2*pi*(0:1/ndir:1-1/ndir)';
dir = [cos(phi),sin(phi)];

qr1 = gauleg(0,1,12);
% qr2 = Duffy(TProd(gauleg(0,1,12)));

% Construct mesh
Mesh.Coordinates = [-1 -1; 1 -1; 1 1; -1 1];
% Mesh.Elements = [1 2 3; 1 4 3];
Mesh.Elements = [1 2 3 4];
Mesh = add_Edges(Mesh);         
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
Mesh.BdFlags(Loc) = -1;
for j = 1:nref
  Mesh = refine_REG(Mesh);
end
Mesh = orient_Elems(Mesh);
Mesh = add_Edge2Elem(Mesh);
Mesh = add_DGData(Mesh);

% Flux parameters
flux_params = {'a',@(omega,h,varargin) 2/(h*omega),...
    'b',@(omega,h,varargin) 0.01*h*omega};
  
% Set PWDG data
Mesh = set_Data_PWDG(Mesh,'Dir',dir,'Omega',omega,flux_params{:});

% Assemble stiffness matrix
[I_Inn,J_Inn,A_Inn] = assemMat_Inn_PDG2(Mesh,@STIMA_Helm_Inn_PWDG,omega);
[I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2(Mesh,[],@STIMA_Helm_Imp_Bnd_PWDG,omega);
A = sparse([I_Inn;I_Bnd],[J_Inn;J_Bnd],[A_Inn;A_Bnd]);

% Assemble load vector
b = assemLoad_Bnd_PDG2(Mesh,[],@LOAD_Imp_Bnd_PWDG,qr1,omega,g);

% Solve equation
u = A\b;

% % Plot solution
% plot_PWDG(u,Mesh,omega,max(0,5-nref));
% return

% Construct evaluation points
nx = 200;
x = [linspace(-1,1,nx)',repmat(0,nx,1)];

% Evaluate numerical solution
u_val = eval_PWDG(Mesh,omega,u,x);

% Evaluate exact solution
u_ex_val = u_ex(x);

% Determine wave number of modified equation
% omega_modeq = dispinv(ndir,omega*2^-nref*d,params_in,3)*2^nref;
% if(real(omega_modeq) < 0)
%   omega_modeq = -omega_modeq;
% end
% omega_modeq = 2*omega-omega_modeq;

% omega_modeq = dispersion(ndir,omega*2^-nref,pi,params_in,3)*2^nref;
% c = eval_PWDG(Mesh,omega,dir,u,[0 0]);
% 
% u_modeq_val = c*exp(i*(omega_modeq*x*d'));

% Determine bloch wave in kernel
[omega_bw,sv_bw,u_loc_bw] = dispersion(ndir,omega*2^-nref,pi,flux_params,4);
omega_bw = omega_bw*2^nref;
x0 = Mesh.Coordinates(Mesh.Elements(:,1),:)*d';
expx0 = exp(i*omega_bw*x0);
expx0 = repmat(expx0(:),1,ndir).';
u_bw = expx0(:).*repmat(u_loc_bw,size(Mesh.Elements,1),1);
u_bw_val = eval_PWDG(Mesh,omega,u_bw,x);
u_bw_val_0 = eval_PWDG(Mesh,omega,u_bw,0);
u_bw_val = u_bw_val/(max(abs(u_bw_val))*u_bw_val_0/abs(u_bw_val_0))*(u_ex([0 0])/abs(u_ex([0 0])));

% % Plot Bloch wave
% plot_PWDG(u_bw,Mesh,omega,max(0,5-nref));

% Plot solutions
figure;
plot(x(:,1),real(u_val),'-',x(:,1),real(u_ex_val),'-',x(:,1),real(u_bw_val),'-');
hold on;
plot(x(:,1),imag(u_val),':',x(:,1),imag(u_ex_val),':',x(:,1),imag(u_bw_val),':');
hold off;
legend('numerical','exact','bloch wave');
xlabel('\bf x');
ylabel('\bf u');
title('\bf Numerical PWDG and exact solutions of 2D Helmholtz equation');

return