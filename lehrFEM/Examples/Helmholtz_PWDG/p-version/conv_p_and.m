function [] = conv_p_and()
% conv_p_and.m                                  5.5.2009
% p-convergence for the PWDG, different domains, meshes, plots, norms, solutions  
% It should work to check the order of convergence for singular and regular
% solutions
%
% plot errors:
% L^2-proj, UWVFsolution, PWDGsolution, in  _ L^2 norm
%                                           _ broken H^1 seminorm
%                                           _ L^2 norm of jumps
% in _ semilogy plot wrt P
%    _ loglog plot wrt (P/ log P)
%
% Finally it plot the mesh
%
% (results exponential for regular solutions, strange for singular,
% algebraic?)
%
% It need MultiNormErr_PWDG.m and ...
%
%
% Andrea Moiola, modified from conv_p_l2_0.m
%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

tic
% Define problem parameters
BC=2;       % =1 DIRICHLET BC, =2 IMPEDENCE BC

omega = 10                     % omega = wavenumber of approximating plane waves
d = [1 0.5];d = d/norm(d);
lambda = omega*1;               %lambda = wavenumber in u_ex, = omega for Helmholtz solutions, for other cases, set function 'f' below (rhs)
psi = 0;

vnorm=@(x) sqrt(x(:,1).^2+x(:,2).^2) ;          %euclid norm 2D for arrays N x 2 (--> result Nx1)
vangle = @(x) angle(x(:,1)+i*x(:,2));

%% exact solutions
% PLANE WAVE
u_ex = @(x,varargin) sin(lambda*x*d'+psi);
% % gI = @(x,n,varargin) lambda*(d*n')*cos(lambda*x*d'+psi) + i*omega*sin(lambda*x*d'+psi);
  
Dx_u_ex= @(x,varargin) lambda * d(1)* cos(lambda*x*d'+psi);
Dy_u_ex= @(x,varargin) lambda * d(2)* cos(lambda*x*d'+psi);


% V_0
%   u_ex=@(x,varargin)(besselj(0, lambda*vnorm(x) )   );
%   gI=@(x,n,varargin) (-lambda*besselj(1 , lambda*vnorm(x) ) .*(x*n') )  ./vnorm(x)  +i*omega*besselj(0, lambda*vnorm(x) );

% V(z^nn) normalized
% nn=1;
% ii=sqrt(-1);
% u_ex = @(x,varargin) exp(ii*nn*angle(x(1)+ii*x(2))) * besselj(nn, lambda*norm(x));
% gI = @(x,n,varargin) (lambda* exp(ii*nn*angle(x(1)+ii*x(2))) * (besselj(nn-1,(lambda*norm(x))) - besselj(nn+1,(lambda*norm(x)))) * x/(2*norm(x)) + ...
%     ii*nn/norm(x)^2 * exp(ii*nn*angle(x(:,1)+ii*x(:,2))) * besselj(nn, lambda*norm(x)) * [-x(2), x(1)] ) * n' + ...
%     ii * omega * exp(ii*nn*angle(x(1)+ii*x(2))) * besselj(nn, lambda*norm(x));


% %sum of 2 PW, ok
% d1 = [1 0.5];   d1 = d1/norm(d1);   d2 = [1 -0.75]; d2 = d2/norm(d2); psi1 = 0;   psi2 = 2;
%   u_ex = @(x,varargin) sin(lambda*x*d1'+psi1) + sin(lambda*x*d2'+psi2);
%   gI = @(x,n,varargin) lambda*(d1*n')*cos(lambda*x*d1'+psi1) + i*omega*sin(lambda*x*d1'+psi1) + lambda*(d2*n')*cos(lambda*x*d2'+psi2) + i*omega*sin(lambda*x*d2'+psi2);

% % hankel, like bessel
% u_ex = @(x,varargin) besselj(0, lambda*vnorm([x(:,1)+1.25,x(:,2)]) )+ sqrt(-1)*bessely(0, lambda*vnorm([x(:,1)+1.25,x(:,2)]));
 

%circular wave for  xi integer,
%Helmholtz corner singularity in the origin for real xi non integer
 %theta_0=pi/4;                      %theta = pi/3 pi/4 for Dirichlet homogeneous solutions on appropriate meshes
 %theta_0=pi/3;
 theta_0=0;
 xi=1                                                                      %xi= 1, 2/3, 3/2
 u_ex = @(x,varargin) (besselj(xi, lambda*vnorm(x) )) .*cos(xi*(vangle(x)-theta_0));
  
 Dx_u_ex= @(x,varargin) cos(vangle(x)) .* cos(xi*(vangle(x)-theta_0))*lambda/2 .* (besselj(xi-1, lambda*vnorm(x)) - besselj(xi+1, lambda*vnorm(x) ) ) + ...
     xi* sin(vangle(x)) .* sin(xi*(vangle(x)-theta_0)) ./ vnorm(x) .* besselj(xi,  lambda*vnorm(x));
 Dy_u_ex= @(x,varargin) sin(vangle(x)) .* cos(xi*(vangle(x)-theta_0))*lambda/2 .* (besselj(xi-1, lambda*vnorm(x)) - besselj(xi+1, lambda*vnorm(x) ) ) - ...
     xi* cos(vangle(x)) .* sin(xi*(vangle(x)-theta_0)) ./ vnorm(x) .* besselj(xi,  lambda*vnorm(x));


% impedence and dirichlet condition if not already defined
gI = @(x,n,varargin) [Dx_u_ex(x), Dy_u_ex(x)] * n'  + i*omega* u_ex(x);
gD=u_ex;


%gD = @(x,varargin) sin(lambda*x*d'+psi);
%f = @(x,varargin) (lambda^2-omega^2)*sin(lambda*x*d'+psi);         %source term when using PW with wave# lambda != omega, Helmholtz non-homogeneous
f= @(x,varargin) 0;

grad_u_ex= @(x,varargin) [Dx_u_ex(x,varargin), Dy_u_ex(x,varargin) ];

%% Define discretization parameters                 +++++++++++++++++++++++

nref = 1;                             % # mesh refinements = 1,2,3 for rectangle, 0 for L-shaped (built already refined)
ndir = 3:27;                        %# PW directions used (=p)
num = numel(ndir);

% Define flux parameters
flcoef=10          %coefficient to improve the convergence of the method with PWDG fluxes
flux_params_pwdg = {'a',@(omega,h,p,varargin) flcoef*p/(h*omega*log(p)),   'b',@(omega,h,p,varargin) h*omega*log(p)/(p*flcoef),'d',@(omega,h,p,varargin) h*omega*log(p)/(p*flcoef)};
flux_params_uwvf = {'a',0.5,'b',0.5,'d',0.5};

% Define quadrature rules (then I can decide where to use qr2 and qr3)
qr1 = gauleg(0,1,12);                   % 1D for edges
qr2 = Duffy(TProd(gauleg(0,1,12)));     % for triangular elements
qr3 = Duffy(TProd(gauleg(0,1,50)));     % for triangular elements  high order to resolve the singularities
% qr2 = TProd(gauleg(0,1,12)); % for quadrilateral elements



%% Construct mesh
% square / rectangle
R=1
offsetR=0;                                                                                               % offset to move the mesh in order to have the singularity on a edge, not on a corner
Mesh.Coordinates = [0 -R/2 + offsetR; R -R/2 + offsetR; R R/2 + offsetR; 0 R/2 + offsetR];               %right square
% Mesh.Coordinates = [-R 0; R 0; R R; -R R];               %upper rectangle
% Mesh.Coordinates = [-0.5 -1; 1 -1; 1 1; -0.5 1];         %right rectangle extended to left
Mesh.Elements = [1 2 3; 1 4 3]; % triangular elements
% Mesh.Elements = [1 2 3 4]; % quadrilateral elements



% % Construct mesh, L-shaped domain

% BBOX=[-1 -1; 1 1];
% h0=0.5;                                                                    %meshsize
% DHandle =@(x) dist_diff( dist_rect(x, [-1,-1],2,2)   ,   dist_rect(x, [-1,-1],1,1) );        % L-shaped around 0
% HHandle= @h_uniform;
% disp=1;                                                                     % =1 show mesh generation, 0 = hide
% FixedPos=[0 0; 0 -1; 1 -1; 1 1; -1 1; -1 0]; 
% Mesh=init_Mesh(BBOX, h0, DHandle, HHandle, FixedPos, disp);
% nref=0;


% % % Construct mesh, unstructured on a square domain
% BBOX=[0 -1; 2 1];
% h0=0.66;                                                                    %meshsize
% DHandle =@(x) dist_rect(x, [0,-1],2,2) ;        % L-shaped around 0
% HHandle= @h_uniform;
% disp=1;                                                                     % =1 show mesh generation, 0 = hide
% FixedPos=[0 0; 0 -1; 2 -1; 2 1; 0 1]; 
% Mesh=init_Mesh(BBOX, h0, DHandle, HHandle, FixedPos, disp);
% nref=0;

% %triangular mesh, corner in  0, 120 deg, built to test homogeneous
% %Dirichlet BC with \xi=3/2
% BBOX=[-1 0; 1 1];
% h0=0.2;                                                                    %meshsize
% DHandle =@(x)  dist_tri(x, [0; 0],[1; 0],[-0.5; sqrt(3)/2]);        
% HHandle= @h_uniform;
% disp=1;                                                                     % =1 show mesh generation, 0 = hide
% FixedPos=[0 0; 0 -1; 1 -1; 1 1; -1 1; -1 0]; 
% Mesh=init_Mesh(BBOX, h0, DHandle, HHandle, FixedPos, disp);
% nref=0;





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
err_pwdg_h1 = nan(1,num);
err_uwvf_h1 = nan(1,num);
err_proj_h1 = nan(1,num);
err_pwdg_dg = nan(1,num);
err_uwvf_dg = nan(1,num);
err_proj_dg = nan(1,num);
%err_proj_M = nan(1,num);    %in L^infty

%% Loop over p
for j=1:num
  
  % Define plane wave directions
  Mesh = set_Data_PWDG(Mesh,'nDofs',ndir(j),'NewDir','replace','Omega',omega);
  
  % Set flux parameters for plane wave DG
  Mesh = set_Data_PWDG(Mesh,'Omega',omega,flux_params_pwdg{:});
  
  % Assemble stiffness matrix for plane waves
  [I_Inn,J_Inn,A_Inn] = assemMat_Inn_PDG2(Mesh,@STIMA_Helm_Inn_PWDG,omega);
  if BC==2
      [I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2(Mesh,[],@STIMA_Helm_Imp_Bnd_PWDG,omega);
  else
      [I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2(Mesh,[],@STIMA_Helm_Dir_Bnd_PWDG,omega);       
  end
  A_pwdg = sparse([I_Inn;I_Bnd],[J_Inn;J_Bnd],[A_Inn;A_Bnd]);

  % Assemble load vector for plane waves
  if BC==2
    b_g = assemLoad_Bnd_PDG2(Mesh,[],@LOAD_Imp_Bnd_PWDG,qr1,omega,gI);
  else
      b_g = assemLoad_Bnd_PDG2(Mesh,[],@LOAD_Dir_Bnd_PWDG,qr1,omega,gD);  
  end
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
  if BC==2  
    [I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2(Mesh,[],@STIMA_Helm_Imp_Bnd_PWDG,omega);
  else
    [I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2(Mesh,[],@STIMA_Helm_Dir_Bnd_PWDG,omega);
  end
  A_uwvf = sparse([I_Inn;I_Bnd],[J_Inn;J_Bnd],[A_Inn;A_Bnd]);

  % Assemble load vector for ultraweak variational formulation
  if BC==2  
    b_g = assemLoad_Bnd_PDG2(Mesh,[],@LOAD_Imp_Bnd_PWDG,qr1,omega,gI);
  else
    b_g = assemLoad_Bnd_PDG2(Mesh,[],@LOAD_Dir_Bnd_PWDG,qr1,omega,gD);
  end
  b_f = assemLoad_Vol_PDG2(Mesh,@LOAD_Vol_PWDG,qr2,omega,f);
  b_uwvf = b_f + b_g;
  
  % Solve equation for ultraweak variational formulation
%   u_uwvf = A_uwvf\b_uwvf;
  [PL,PR] = assemPrec_SVD_PDG(Mesh,A_uwvf);
  w = (PL*A_uwvf*PR)\(PL*b_uwvf);
  u_uwvf = PR*w;
  
  
  
  
  
  % Assemble mass matrix for L^2-projection
  [I_mass_inn,J_mass_inn,M_inn] = assemMat_Inn_PDG2(Mesh,@MASS_Inn_PWDG,omega);
  [I_mass_bnd,J_mass_bnd,M_bnd] = assemMat_Bnd_PDG2(Mesh,[],@MASS_Bnd_PWDG,omega);
  I = [I_mass_inn;I_mass_bnd];
  J = [J_mass_inn;J_mass_bnd];
  M = sparse(I,J,[M_inn;M_bnd]);
  
  % Calculate L2 projection of exact solution onto discretization space
  b_proj = assemLoad_Vol_PDG2(Mesh,@LOAD_Vol_PWDG,qr2,omega,u_ex);
  u_proj = M\b_proj;
    
 
  [err_pwdg(j), err_pwdg_h1(j), err_pwdg_dg(j)] = MultiNormErr_PWDG(Mesh,u_pwdg,qr2,omega,u_ex, grad_u_ex    ,1);           %quadratures! ...
  [err_uwvf(j), err_uwvf_h1(j), err_uwvf_dg(j)] = MultiNormErr_PWDG(Mesh,u_uwvf,qr2,omega,u_ex, grad_u_ex    ,1);
  [err_proj(j), err_proj_h1(j), err_proj_dg(j)] = MultiNormErr_PWDG(Mesh,u_proj,qr2,omega,u_ex, grad_u_ex    ,1);

  
end

err_pwdg= sqrt(err_pwdg);
err_pwdg_h1= sqrt(err_pwdg_h1);
err_pwdg_dg = sqrt(err_pwdg_dg);
err_uwvf = sqrt(err_uwvf);
err_uwvf_h1 = sqrt(err_uwvf_h1);
err_uwvf_dg = sqrt(err_uwvf_dg);
err_proj = sqrt(err_proj);
err_proj_h1 = sqrt(err_proj_h1);
err_proj_dg = sqrt(err_proj_dg);


ndir_l=ndir./log(ndir);             % scaled number of plane waves, to see the order of convergence
% Plot errors
% L2 - semilog
figure;
set(gca, 'FontSize', 18);
plot(ndir,err_pwdg,'-.',ndir,err_uwvf,'--',ndir,err_proj,'-','linewidth',2);
set(gca,'YScale','log','XLim',[min(ndir),max(ndir)],'XTick',ndir(1:2:end));
grid on;
xlabel('\bf number of local plane wave basis functions');
ylabel('\bf L^2 error');
%title(sprintf('%s Convergence of PWDG, h = %g, xi = %d, omega = %d','\bf',2^(1-nref), xi, omega));
legend('PWDG L2','ultra-weak L2 ','proj. L2');

% H1 - semilog
figure;
set(gca, 'FontSize', 18);
plot(ndir,err_pwdg_h1,'-.',ndir,err_uwvf_h1,'--',ndir,err_proj_h1,'-','linewidth',2);
set(gca,'YScale','log','XLim',[min(ndir),max(ndir)],'XTick',ndir(1:2:end));
grid on;
xlabel('\bf number of local plane wave basis functions');
ylabel('\bf broken H^1 error');
%title(sprintf('%s Convergence of PWDG for ex. sol. in kernel, h = %g, xi = %d, omega = %d','\bf',2^(1-nref), xi, omega));
legend('PWDG br. H^1','ultra-weak br. H^1.','proj. br. H^1');

% Jumps - semilog
figure;
set(gca, 'FontSize', 18);
plot(ndir,err_pwdg_dg,'-.',ndir,err_uwvf_dg,'--',ndir,err_proj_dg,'-','linewidth',2);
set(gca,'YScale','log','XLim',[min(ndir),max(ndir)],'XTick',ndir(1:2:end));
grid on;
xlabel('\bf number of local plane wave basis functions');
ylabel('\bf jumps L^2-norm');
%title(sprintf('%s Convergence of PWDG for ex. sol. in kernel, h = %g, xi = %d, omega = %d','\bf',2^(1-nref), xi, omega));
legend('PWDG jumps','ultra-weak jumps ','proj. jumps');



% L2 - loglog
figure;
set(gca, 'FontSize', 18);
loglog(ndir_l,err_pwdg,'-.',ndir_l,err_uwvf,'--',ndir_l,err_proj,'-','linewidth',2);
grid on;
xlabel('\bf p/log(p)');
ylabel('\bf L^2 error');
%title(sprintf('%s Convergence of PWDG, h = %g, xi = %d, omega = %d','\bf',2^(1-nref), xi, omega));
legend('PWDG L2','ultra-weak L2','proj. L2');

% H1 - loglog
figure;
set(gca, 'FontSize', 18);
loglog(ndir_l,err_pwdg_h1,'-.',ndir_l,err_uwvf_h1,'--',ndir_l,err_proj_h1,'-','linewidth',2);
%set(gca,'YScale','log','XLim',[min(ndir),max(ndir)],'XTick',ndir);
grid on;
xlabel('\bf p/log(p)');
ylabel('\bf broken H^1 error');
%title(sprintf('%s Convergence of PWDG for ex. sol. in kernel, h = %g, xi = %d, omega = %d','\bf',2^(1-nref), xi, omega));
legend('PWDG br. H^1','ultra-weak br. H^1.','proj. br. H^1');

% Jumps - loglog
figure;
set(gca, 'FontSize', 18);
loglog(ndir_l,err_pwdg_dg,'-.',ndir_l,err_uwvf_dg,'--',ndir_l,err_proj_dg,'-','linewidth',2);
%set(gca,'YScale','log','XLim',[min(ndir),max(ndir)],'XTick',ndir);
grid on;
xlabel('\bf p/log(p)');
ylabel('\bf jumps L^2-norm');
%title(sprintf('%s Convergence of PWDG for ex. sol. in kernel, h = %g, xi = %d, omega = %d','\bf',2^(1-nref), xi, omega));
legend('PWDG jumps','ultra-weak jumps','proj. jumps');


plot_Mesh(Mesh);  set(gca, 'FontSize', 18); axis on

toc


%add_Slope(gca,'SouthWest',-5.5);

% for the regular function in paper, pictures adjusted with something like:
% figure(1); axis([3 27 10^(-12) 1]);   set(gca, 'ytick', 10.^[-12:2:0]);
% figure(2); axis([3 27 10^(-10) 100]); set(gca, 'ytick', 10.^[-10:2:2]);

