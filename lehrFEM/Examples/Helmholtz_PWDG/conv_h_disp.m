function [] = conv_h_disp()
%CONV_H_DISP compare convergence for more and less dispersive directions
%   Detailed explanation goes here

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Define nondispersive problem
omega = 10;
d = [1 0.1];
d = d/norm(d);
lambda = omega;
psi = 0;
u_ex_0 = @(x,varargin) sin(lambda*x*d'+psi);
du_ex_0 = @(x,varargin) lambda*cos(lambda*x*d([1 1],:)'+psi).*d(ones(size(x,1),1),:);
g_0 = @(x,n,varargin) lambda*(d*n')*cos(lambda*x*d'+psi) + i*omega*sin(lambda*x*d'+psi);
f_0 = @(x,varargin) (lambda^2-omega^2)*sin(lambda*x*d'+psi);

% Define dispersive problem (for 4 plane wave shape functions)
d = [1 1];
d = d/norm(d);
u_ex_1 = @(x,varargin) sin(lambda*x*d'+psi);
du_ex_1 = @(x,varargin) lambda*cos(lambda*x*d([1 1],:)'+psi).*d(ones(size(x,1),1),:);
g_1 = @(x,n,varargin) lambda*(d*n')*cos(lambda*x*d'+psi) + i*omega*sin(lambda*x*d'+psi);
f_1 = @(x,varargin) (lambda^2-omega^2)*sin(lambda*x*d'+psi);

% Define discretization parameters
nref = [1 5];
num = nref(2)-nref(1)+1;
ndir = 4;
phi = 2*pi*(0:1/ndir:1-1/ndir)';
dir = [cos(phi),sin(phi)];

% Define flux parameters
flux_params = {'a',@(omega,h,varargin)1/(h*omega),...
  'b',@(omega,h,varargin)0.1*h*omega,'d',1};

% Define quadrature rules
qr1 = gauleg(0,1,12);
qr2 = Duffy(TProd(gauleg(0,1,12)));

% Calculate mesh width
h = 2.^(1-nref(1):-1:1-nref(2));

% Initialize errors
err_pwdg_0 = nan(1,num);
err_pwdg_1 = nan(1,num);
err_proj_0 = nan(1,num);
err_proj_1 = nan(1,num);

% Initialize mesh
Mesh.Coordinates = [-1 -1; 1 -1; 1 1; -1 1];
Mesh.Elements = [1 2 3; 1 4 3];
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
  Mesh = set_Data_PWDG(Mesh,'Dir',dir,'NewDir','replace','Omega',omega,flux_params{:});
  
  % Assemble stiffness matrix
  [I_Inn,J_Inn,A_Inn] = assemMat_Inn_PDG2(Mesh,@STIMA_Helm_Inn_PWDG,omega);
  [I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2(Mesh,[],@STIMA_Helm_Imp_Bnd_PWDG,omega);
  A_pwdg = sparse([I_Inn;I_Bnd],[J_Inn;J_Bnd],[A_Inn;A_Bnd]);

  % Assemble load vector for nondispersive problem
  b_g = assemLoad_Bnd_PDG2(Mesh,[],@LOAD_Imp_Bnd_PWDG,qr1,omega,g_0);
  b_f = assemLoad_Vol_PDG2(Mesh,@LOAD_Vol_PWDG,qr2,omega,f_0);
  b_pwdg_0 = b_f + b_g;
  
  % Solve equation for nondispersive problem
  u_pwdg_0 = A_pwdg\b_pwdg_0;
  
  % Assemble load vector for dispersive problem
  b_g = assemLoad_Bnd_PDG2(Mesh,[],@LOAD_Imp_Bnd_PWDG,qr1,omega,g_1);
  b_f = assemLoad_Vol_PDG2(Mesh,@LOAD_Vol_PWDG,qr2,omega,f_1);
  b_pwdg_1 = b_f + b_g;
  
  % Solve equation for dispersive problem
  u_pwdg_1 = A_pwdg\b_pwdg_1;
  
  % Assemble energy-norm inner product
  [I_h1_inn,J_h1_inn,H_inn] = assemMat_Inn_PDG2(Mesh,@SCAPRO_Helm_Inn_PWDG,omega);
  [I_h1_bnd,J_h1_bnd,H_bnd] = assemMat_Bnd_PDG2(Mesh,[],@SCAPRO_Helm_Bnd_PWDG,omega);
  I_h1 = [I_h1_inn;I_h1_bnd];
  J_h1 = [J_h1_inn;J_h1_bnd];
  H_h1 = [H_inn;H_bnd];
  B = sparse(I_h1,J_h1,H_h1);
  
  % Calculate projection for nondispersive problem
  b_proj_0 = assemLoad_Vol_PDG2(Mesh,@LOAD_EnergyProj_Vol_PWDG,qr2,omega,u_ex_0,du_ex_0);
  u_proj_0 = B\b_proj_0;
  
  % Calculate projection for dispersive problem
  b_proj_1 = assemLoad_Vol_PDG2(Mesh,@LOAD_EnergyProj_Vol_PWDG,qr2,omega,u_ex_1,du_ex_1);
  u_proj_1 = B\b_proj_1;
  
  % Calculate errors
  err_pwdg_0(j) = EnergyErr_PWDG(Mesh,u_pwdg_0,qr2,omega,u_ex_0,du_ex_0,0);
  err_pwdg_1(j) = EnergyErr_PWDG(Mesh,u_pwdg_1,qr2,omega,u_ex_1,du_ex_1,0);
  err_proj_0(j) = EnergyErr_PWDG(Mesh,u_proj_0,qr2,omega,u_ex_0,du_ex_0,0);
  err_proj_1(j) = EnergyErr_PWDG(Mesh,u_proj_1,qr2,omega,u_ex_1,du_ex_1,0);
  
end

% Plot errors
figure;
plot(h,err_pwdg_0,'--',h,err_pwdg_1,'-');
hold on;
plot(h,err_proj_0,':',h,err_proj_1,':');
hold off;
set(gca,'XScale','log','XDir','reverse','YScale','log','XLim',[min(h),max(h)]);
grid on;
xlabel('\bf h');
ylabel('\bf energy-norm error');
title(sprintf('%s Effect of dispersion on convergence','\bf'));
legend('sol. for nondisp. prob.','sol. for disp. prob.','proj. for nondisp. prob.','proj. for disp. prob.','Location','SouthWest');

return