function [] = conv_p_l2_1()
%CONV_P_L2_1 p-convergence

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Define problem parameters
omega = 1;
d = [1 0.5];
d = d/norm(d);
lambda = 0.5*omega;
psi = 0;
u_ex = @(x,varargin) sin(lambda*x*d'+psi);

gI = @(x,n,varargin) lambda*(d*n')*cos(lambda*x*d'+psi) + i*omega*sin(lambda*x*d'+psi);
% gD = @(x,varargin) sin(lambda*x*d'+psi);
f = @(x,varargin) (lambda^2-omega^2)*sin(lambda*x*d'+psi);

% Define discretization parameters
nref = 2;
ndir = 3:12;
num = numel(ndir);

% Define flux parameters
flux_params_pwdg = {'a',@(omega,h,varargin) 1/(h*omega),'b',0,'d',1};
% flux_params_pwdg = {'a',@(omega,h,varargin) 1/(h*omega),...
%   'b',@(omega,h,varargin) 0.01*h*omega,'d',0.5};
flux_params_uwvf = {'a',0.5,'b',0.5,'d',0.5};

% Define quadrature rules
qr1 = gauleg(0,1,12);
qr2 = Duffy(TProd(gauleg(0,1,12))); % for triangular elements
% qr2 = TProd(gauleg(0,1,12)); % for quadrilateral elements

% Construct mesh
Mesh.Coordinates = [-1 -1; 1 -1; 1 1; -1 1];
Mesh.Elements = [1 2 3; 1 4 3]; % triangular elements
% Mesh.Elements = [1 2 3 4]; % quadrilateral elements

Mesh = add_Edges(Mesh);         
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
Mesh.BdFlags(Loc) = -1;
for j=1:nref
  Mesh = refine_REG(Mesh);
end
Mesh = orient_Elems(Mesh);
Mesh = add_Edge2Elem(Mesh);
Mesh = add_DGData(Mesh);

% Initialize errors
err_pwdg = nan(1,num);
err_uwvf = nan(1,num);
err_proj = nan(1,num);

% Loop over meshes
for j=1:num
  
  % Define plane wave directions
  Mesh = set_Data_PWDG(Mesh,'nDofs',ndir(j),'NewDir','replace','Omega',omega);
  
  % Set flux parameters for plane wave DG
  Mesh = set_Data_PWDG(Mesh,'Omega',omega,flux_params_pwdg{:});
  
  % Assemble stiffness matrix for plane waves
  [I_Inn,J_Inn,A_Inn] = assemMat_Inn_PDG2(Mesh,@STIMA_Helm_Inn_PWDG,omega);
  [I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2(Mesh,[],@STIMA_Helm_Imp_Bnd_PWDG,omega);
%   [I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2(Mesh,[],@STIMA_Helm_Dir_Bnd_PWDG,omega);
  A_pwdg = sparse([I_Inn;I_Bnd],[J_Inn;J_Bnd],[A_Inn;A_Bnd]);

  % Assemble load vector for plane waves
  b_g = assemLoad_Bnd_PDG2(Mesh,[],@LOAD_Imp_Bnd_PWDG,qr1,omega,gI);
%   b_g = assemLoad_Bnd_PDG2(Mesh,[],@LOAD_Dir_Bnd_PWDG,qr1,omega,gD);
  b_f = assemLoad_Vol_PDG2(Mesh,@LOAD_Vol_PWDG,qr2,omega,f);
  b_pwdg = b_f + b_g;
  
  % Solve equation for plane waves
%   u_pwdg = A_pwdg\b_pwdg;
  [PL,PR] = assemPrec_SVD_PDG(Mesh,A_pwdg);
  w = (PL*A_pwdg*PR)\(PL*b_pwdg);
  u_pwdg = PR*w;
  
  % Set flux parameters for ultraweak variational formulation
  Mesh = set_Data_PWDG(Mesh,'Omega',omega,flux_params_uwvf{:});
  
  % Assemble stiffness matrix for ultraweak variational formulation
  [I_Inn,J_Inn,A_Inn] = assemMat_Inn_PDG2(Mesh,@STIMA_Helm_Inn_PWDG,omega);
  [I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2(Mesh,[],@STIMA_Helm_Imp_Bnd_PWDG,omega);
%   [I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2(Mesh,[],@STIMA_Helm_Dir_Bnd_PWDG,omega);
  A_uwvf = sparse([I_Inn;I_Bnd],[J_Inn;J_Bnd],[A_Inn;A_Bnd]);

  % Assemble load vector for ultraweak variational formulation
  b_g = assemLoad_Bnd_PDG2(Mesh,[],@LOAD_Imp_Bnd_PWDG,qr1,omega,gI);
%   b_g = assemLoad_Bnd_PDG2(Mesh,[],@LOAD_Dir_Bnd_PWDG,qr1,omega,gD);
  b_f = assemLoad_Vol_PDG2(Mesh,@LOAD_Vol_PWDG,qr2,omega,f);
  b_uwvf = b_f + b_g;
  
  % Solve equation for ultraweak variational formulation
%   u_uwvf = A_uwvf\b_uwvf;
  [PL,PR] = assemPrec_SVD_PDG(Mesh,A_uwvf);
  w = (PL*A_uwvf*PR)\(PL*b_uwvf);
  u_uwvf = PR*w;
  
  % Assemble mass matrix
  [I_mass_inn,J_mass_inn,M_inn] = assemMat_Inn_PDG2(Mesh,@MASS_Inn_PWDG,omega);
  [I_mass_bnd,J_mass_bnd,M_bnd] = assemMat_Bnd_PDG2(Mesh,[],@MASS_Bnd_PWDG,omega);
  I = [I_mass_inn;I_mass_bnd];
  J = [J_mass_inn;J_mass_bnd];
  M = sparse(I,J,[M_inn;M_bnd]);
  
  % Calculate L2 projection of exact solution onto discretization space
  b_proj = assemLoad_Vol_PDG2(Mesh,@LOAD_Vol_PWDG,qr2,omega,u_ex);
  u_proj = M\b_proj;
  
  % Calculate errors
  err_pwdg(j) = L2Err_PWDG(Mesh,u_pwdg,qr2,omega,u_ex);
  err_uwvf(j) = L2Err_PWDG(Mesh,u_uwvf,qr2,omega,u_ex);
  err_proj(j) = L2Err_PWDG(Mesh,u_proj,qr2,omega,u_ex);
  
end

% Plot errors
figure;
plot(ndir,err_pwdg,'-',ndir,err_uwvf,'--',ndir,err_proj,':');
set(gca,'YScale','log','XLim',[min(ndir),max(ndir)],'XTick',ndir);
grid on;
xlabel('\bf number of local plane wave basis functions');
ylabel('\bf L^2 error');
title(sprintf('%s Convergence of PWDG, h = %g','\bf',2^(1-nref)));
legend('PWDG','ultra-weak','proj.');
