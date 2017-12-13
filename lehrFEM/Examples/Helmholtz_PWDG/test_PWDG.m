function [] = test_PWDG()
%TEST_PWDG solve sample problem with PWDG

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% should print
%   L2 Error : 0.00662562    L2 Proj.Error : 0.00299635
%   En.Error : 0.365672    En.Proj.Error : 0.35686

% Define problem parameters
omega = 1.5*pi;
d = [1 0.5];
d = d/norm(d);
lambda = 0.5*omega;
psi = 0;
u_ex = @(x,varargin) sin(lambda*x*d'+psi);
du_ex = @(x,varargin) lambda*cos(lambda*x*d([1 1],:)'+psi).*d(ones(size(x,1),1),:);
gI = @(x,n,varargin) lambda*(d*n')*cos(lambda*x*d'+psi) + i*omega*sin(lambda*x*d'+psi);
% gD = @(x,varargin) sin(lambda*x*d'+psi);
f = @(x,varargin) (lambda^2-omega^2)*sin(lambda*x*d'+psi);

% omega = 0.5;
% d = [5/4 -i*3/4];
% u_ex = @(x,varargin) exp(i*omega*x*d.');
% du_ex = @(x,varargin) i*omega*exp(i*omega*x*d([1 1],:).').*d(ones(size(x,1),1),:);
% gI = @(x,n,varargin) i*omega*(1+d*n.')*exp(i*omega*x*d.');
% f = @(x,varargin) zeros(size(x,1),1);

% Define discretization parameters
nref = 4;
ndir = 7;
phi = 2*pi*(0:1/ndir:1-1/ndir)';
dir = [cos(phi),sin(phi)];
% ndir = 3;
% phi = 2*pi*(0:1/ndir:1-1/ndir)';
% dir = [cos(phi),sin(phi)];
% dir(1,:) = d;

qr1 = gauleg(0,1,12);
qr2 = Duffy(TProd(gauleg(0,1,12)));

% Construct mesh
Mesh.Coordinates = [-1 -1; 1 -1; 1 1; -1 1];
Mesh.Elements = [1 2 3; 1 4 3];
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
% Mesh = set_Data_PWDG(Mesh,'Dir',dir,'Omega',omega,...
%   'a',@(w,h,varargin) 1/(w*h),'b',0,'d',1);
% Mesh = set_Data_PWDG(Mesh,'Dir',dir,'Omega',omega,...
%   'a',@(w,h,varargin) 2/(w*h),'b',@(w,h,varargin) 0.01*(w*h),'d',0.5);
Mesh = set_Data_PWDG(Mesh,'nDofs',ndir,'MakeDir','random','Omega',omega,...
  'a',@(w,h,varargin) 2/(w*h),'b',@(w,h,varargin) 0.01*w*h);
Mesh = set_Data_PWDG(Mesh,'RemDir',1:ndir);
Mesh = set_Data_PWDG(Mesh,'Omega',omega,'a',@(w,h,varargin) 1/(w*h),'b',0,'d',0);
for j=1:ndir
  Mesh = set_Data_PWDG(Mesh,'Dir',dir(j,:),'DofData',j);
end
% nElem = size(Mesh.Elements,1);
% Mesh = set_Data_PWDG(Mesh,'nDofs',ndir,'MakeDir','orient','Phi0',2*pi*(0:1/nElem:1-1/nElem),'NewDir','replace');

% Assemble stiffness matrix
[I_Inn,J_Inn,A_Inn] = assemMat_Inn_PDG2(Mesh,@STIMA_Helm_Inn_PWDG,omega);
[I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2(Mesh,[],@STIMA_Helm_Imp_Bnd_PWDG,omega);
% [I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2(Mesh,[],@STIMA_Helm_Dir_Bnd_PWDG,omega);
A = sparse([I_Inn;I_Bnd],[J_Inn;J_Bnd],[A_Inn;A_Bnd]);

% Assemble load vector
b_g = assemLoad_Bnd_PDG2(Mesh,[],@LOAD_Imp_Bnd_PWDG,qr1,omega,gI);
% b_g = assemLoad_Bnd_PDG2(Mesh,[],@LOAD_Dir_Bnd_PWDG,qr1,omega,gD);
b_f = assemLoad_Vol_PDG2(Mesh,@LOAD_Vol_PWDG,qr2,omega,f);
b = b_f + b_g;

% % Solve equation
% u = A\b;

% Construct preconditioner
[PL,PR] = assemPrec_SVD_PDG([Mesh.ElemData.nDofs],A);

% Solve equation with preconditioning
w = (PL*A*PR)\(PL*b);
u = PR*w;

% % Solve regularized equation
% t = 1e-5;
% u = (A'*A+t*speye(size(A,1)))\(A'*b);

% % Get l2 norm
% fprintf('l2 norm of solution: %g\n',norm(u));

% % Plot solution
% plot_PWDG(u,Mesh,omega,max(0,5-nref));
% plot_PWDG_Dir(Mesh,u);
% contour_PWDG(u,Mesh);
% plot_Mesh(Mesh,'tas');

% Assemble L2 inner product
[I_mass_inn,J_mass_inn,M_inn] = assemMat_Inn_PDG2(Mesh,@MASS_Inn_PWDG,omega);
[I_mass_bnd,J_mass_bnd,M_bnd] = assemMat_Bnd_PDG2(Mesh,[],@MASS_Bnd_PWDG,omega);
M = sparse([I_mass_inn;I_mass_bnd],[J_mass_inn;J_mass_bnd],[M_inn;M_bnd]);

% Calculate L2 projection of exact solution onto discretization space
b_u = assemLoad_Vol_PDG2(Mesh,@LOAD_Vol_PWDG,qr2,omega,u_ex);
u_p = M\b_u;
% plot_PWDG(u_p,Mesh,omega,max(0,5-nref));

% Calculate L2 discretization error
err = L2Err_PWDG(Mesh,u,qr2,omega,u_ex);
err_p = L2Err_PWDG(Mesh,u_p,qr2,omega,u_ex);
fprintf('L2 Error : %g    L2 Proj.Error : %g\n',err,err_p);

% Assemble energy-norm inner product
[I_h1_inn,J_h1_inn,H_inn] = assemMat_Inn_PDG2(Mesh,@SCAPRO_Helm_Inn_PWDG,omega);
[I_h1_bnd,J_h1_bnd,H_bnd] = assemMat_Bnd_PDG2(Mesh,[],@SCAPRO_Helm_Bnd_PWDG,omega);
B = sparse([I_h1_inn;I_h1_bnd],[J_h1_inn;J_h1_bnd],[H_inn;H_bnd]);

% Compute energy-norm projection of exact solution onto PWDG-space
b_u = assemLoad_Vol_PDG2(Mesh,@LOAD_EnergyProj_Vol_PWDG,qr2,omega,u_ex,du_ex);
u_p = B\b_u;
% plot_PWDG(u_p,Mesh,omega,max(0,5-nref));

% Calculate energy norm discretization error
err = EnergyErr_PWDG(Mesh,u,qr2,omega,u_ex,du_ex,0);
err_p = EnergyErr_PWDG(Mesh,u_p,qr2,omega,u_ex,du_ex,0);
fprintf('En.Error : %g    En.Proj.Error : %g\n\n',err,err_p);

% % Test local changes in basis functions
% Els = [1:35,64,73,47,101];
% dirs = repmat([-1 0],length(Els),1);
% [Mesh1,Data] = set_Data_PWDG(Mesh,'Dir',dirs,'Elements',Els);
% u1 = zeros(size(Mesh1.Elements,1),1);
% u1(Data.I1) = u(Data.I0);
% % plot_PWDG(u1,Mesh1,omega,max(0,5-nref));
% err10 = L2Err_PWDG(Mesh1,u1,qr2,omega,u_ex);
% err1e = EnergyErr_PWDG(Mesh1,u1,qr2,omega,u_ex,du_ex,0);
% fprintf('L2 Error : %g    En.Error : %g\n\n',err10,err1e);

return