function [] = conv_h_en_0()
%CONV_H_EN_0 h-convergence for kernel elements

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Define problem parameters
omega = 1;
% d = [1 0.5];
% d = d/norm(d);
% lambda = omega;
% psi = 0;
% u_ex = @(x,varargin) sin(lambda*x*d'+psi);
% du_ex = @(x,varargin) lambda*cos(lambda*x*d([1 1],:)'+psi).*d(ones(size(x,1),1),:);
% gI = @(x,n,varargin) lambda*(d*n')*cos(lambda*x*d'+psi) + i*omega*sin(lambda*x*d'+psi);
% % gD = @(x,varargin) sin(lambda*x*d'+psi);
% f = @(x,varargin) (lambda^2-omega^2)*sin(lambda*x*d'+psi);

d = [5/4 -i*3/4];
u_ex = @(x,varargin) exp(i*omega*x*d.');
du_ex = @(x,varargin) i*omega*exp(i*omega*x*d([1 1],:).').*d(ones(size(x,1),1),:);
gI = @(x,n,varargin) i*omega*(1+d*n.')*exp(i*omega*x*d.');
f = @(x,varargin) zeros(size(x,1),1);

% Define discretization parameters
nref = [1 5];
num = nref(2)-nref(1)+1;
ndir = 5;
phi = 2*pi*(0:1/ndir:1-1/ndir)';
dir = [cos(phi),sin(phi)];

% Define quadrature rules
qr1 = gauleg(0,1,12);
qr2 = Duffy(TProd(gauleg(0,1,12))); % for triangular elements
% qr2 = TProd(gauleg(0,1,12)); % for quadrilateral elements

% Calculate mesh width
h = 2.^(1-nref(1):-1:1-nref(2));

% Initialize errors
err_pwdg = nan(1,num);
err_uwvf = nan(1,num);
err_proj = nan(1,num);

% Initialize mesh
Mesh.Coordinates = [-1 -1; 1 -1; 1 1; -1 1];
Mesh.Elements = [1 2 3; 1 4 3]; % triangular elements
% Mesh.Elements = [1 2 3 4]; % quadrilateral elements

Mesh = add_Edges(Mesh);         
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
Mesh.BdFlags(Loc) = -1;
for j=1:nref(1)-1
  Mesh = refine_REG(Mesh);
end

% Loop over meshes
for j=1:num
  
  % Refine mesh and add fields
  Mesh = refine_REG(Mesh);
  Mesh = orient_Elems(Mesh);
  Mesh = add_Edge2Elem(Mesh);
  Mesh = add_DGData(Mesh);
  Mesh = set_Data_PWDG(Mesh,'Dir',dir,'Omega',omega);
  
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
  u_uwvf = A_uwvf\b_uwvf;
  
  % Set flux parameters
  Mesh = set_Data_PWDG(Mesh,'Omega',omega,'a',@(w,h,varargin) 1/(w*h),'b',0,'d',1);
  
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
  u_pwdg = A_pwdg\b_pwdg;
  
  % Assemble energy-norm inner product
  [I_h1_inn,J_h1_inn,H_inn] = assemMat_Inn_PDG2(Mesh,@SCAPRO_Helm_Inn_PWDG,omega);
  [I_h1_bnd,J_h1_bnd,H_bnd] = assemMat_Bnd_PDG2(Mesh,[],@SCAPRO_Helm_Bnd_PWDG,omega);
  I_h1 = [I_h1_inn;I_h1_bnd];
  J_h1 = [J_h1_inn;J_h1_bnd];
  H_h1 = [H_inn;H_bnd];
  B = sparse(I_h1,J_h1,H_h1);
  
  % Calculate projection of exact solution onto discretization space
  b_proj = assemLoad_Vol_PDG2(Mesh,@LOAD_EnergyProj_Vol_PWDG,qr2,omega,u_ex,du_ex);
  u_proj = B\b_proj;
  
  % Calculate errors
  err_pwdg(j) = EnergyErr_PWDG(Mesh,u_pwdg,qr2,omega,u_ex,du_ex,0);
  err_uwvf(j) = EnergyErr_PWDG(Mesh,u_uwvf,qr2,omega,u_ex,du_ex,0);
  err_proj(j) = EnergyErr_PWDG(Mesh,u_proj,qr2,omega,u_ex,du_ex,0);
    
end

% Plot errors
figure;
plot(h,err_pwdg,'-',h,err_uwvf,'--',h,err_proj,':');
set(gca,'XScale','log','XDir','reverse','YScale','log','XLim',[min(h),max(h)]);
grid on;
xlabel('\bf h');
ylabel('\bf energy-norm error');
title(sprintf('%s Convergence of PWDG for ex. sol. in kernel, %d local basis fn.','\bf',ndir));
legend('PWDG','ultra-weak','proj.');

% Add convergence rate
p = polyfit(log(h),log(err_pwdg),1);
add_Slope(gca,'SouthEast',p(1));
