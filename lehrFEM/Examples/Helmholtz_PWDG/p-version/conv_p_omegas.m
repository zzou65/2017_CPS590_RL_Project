function [] = conv_p_omegas()
% conv_p_omegas.m                                  5.5.2009
% p-convergence for the PWDG with different wavenumber omegas,
% modification of conv_p_and.m
%
% plots errors only in L^2
% but computed also in broken H^1 seminorm and L^2 norm of jumps
% for L^2-proj, UWVFsolution, PWDGsolution
% in _ semilogy plot wrt P
%    _ loglog plot wrt (P/ log P)
%
%
% It need MultiNormErr_PWDG.m
% The pictures need some manipulation to be put in the paper.
%
% Andrea Moiola, modified from conv_p_l2_0.m
%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

tic
% Define problem parameters
BC=2;       % =1 DIRICHLET BC, =2 IMPEDENCE BC


%omegas=logspace(-1,3,5)
%omegas=[0.3,1,5,10,20, 30,50]
omegas=[0.25,1,4,16,64]
n_om=length(omegas);



%% Define discretization parameters                 +++++++++++++++++++++++

nref = 1;                             % # mesh refinements = 1,2,3 for rectangle, 0 for L-shaped (built already refined)
ndir = 3:1:40;                        %# PW directions used (=p)
num = numel(ndir);

% Define flux parameters
flcoef=10;  %coefficient to improve the convergence of the method with PWDG fluxes
flux_params_pwdg = {'a',@(omega,h,p,varargin) flcoef*p/(h*omega*log(p)),   'b',@(omega,h,p,varargin) h*omega*log(p)/(p*flcoef),'d',@(omega,h,p,varargin) h*omega*log(p)/(p*flcoef)};
flux_params_uwvf = {'a',0.5,'b',0.5,'d',0.5};

% Define quadrature rules (then I can decide where to use qr2 and qr3)
qr1 = gauleg(0,1,12);                   % 1D for edges
qr2 = Duffy(TProd(gauleg(0,1,12)));     % for triangular elements
qr3 = Duffy(TProd(gauleg(0,1,50)));     % for triangular elements  high order to resolve the singularities
% qr2 = TProd(gauleg(0,1,12)); % for quadrilateral elements



% Initialize errors
err_pwdg = nan(n_om,num);
err_uwvf = nan(n_om,num);
err_proj = nan(n_om,num);
err_pwdg_h1 = nan(n_om,num);
err_uwvf_h1 = nan(n_om,num);
err_proj_h1 = nan(n_om,num);
err_pwdg_dg = nan(n_om,num);
err_uwvf_dg = nan(n_om,num);
err_proj_dg = nan(n_om,num);
%err_proj_M = nan(n_om,num);    %in L^infty




%% loop over omegas

for j_om=1:n_om
    omega=omegas(j_om);
    
    d = [1 0.5];
    d = d/norm(d);
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

    %circular wave for  xi integer,
    %Helmholtz corner singularity in the origin for real xi non integer
     theta_0=0;
     xi=1                                                                         %xi= 1, 2/3, 3/2
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

   
    
    %% Construct mesh
    Mesh=[];                            %delete old mesh, --> I can build a new scaled one
    % square / rectangle
    R=1
    offsetR=0;                                                                                               % offset to move the mesh in order to have the singularity on a edge, not on a corner
    Mesh.Coordinates = [0 -R/2 + offsetR; R -R/2 + offsetR; R R/2 + offsetR; 0 R/2 + offsetR];               %right square
    % Mesh.Coordinates = [-R 0; R 0; R R; -R R];               %upper rectangle
    % Mesh.Coordinates = [-0.5 -1; 1 -1; 1 1; -0.5 1];         %right rectangle extended to left
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

      [err_pwdg(j_om,j), err_pwdg_h1(j_om,j), err_pwdg_dg(j_om,j)] = MultiNormErr_PWDG(Mesh,u_pwdg,qr2,omega,u_ex, grad_u_ex    ,1);           %quadratures! ...
      [err_uwvf(j_om,j), err_uwvf_h1(j_om,j), err_uwvf_dg(j_om,j)] = MultiNormErr_PWDG(Mesh,u_uwvf,qr2,omega,u_ex, grad_u_ex    ,1);
      [err_proj(j_om,j), err_proj_h1(j_om,j), err_proj_dg(j_om,j)] = MultiNormErr_PWDG(Mesh,u_proj,qr2,omega,u_ex, grad_u_ex    ,1);


    end
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



%% Plot errors
ndir_l=ndir./log(ndir);             % scaled number of plane waves, to see the order of convergence

% L2 - loglog
figure;
%plot(ndir,err_uwvf,'--',ndir,err_pwdg,'-.');
for j=1:n_om
    plot(ndir,err_uwvf(j,:),'-', 'color',[0,0,1-(j-1)/(n_om-1)],'linewidth',2);
    hold on
end
grid off;
xlabel('\bf number of local plane wave basis functions','fontsize',18);
ylabel('\bf L^2 error','fontsize',18);
set(gca,'YScale','log','XLim',[min(ndir),max(ndir)]);% ,'XTick',4:2:ndir(end));
%title(sprintf('%s Convergence of PWDG, h = %g, xi = %d, omega = %d','\bf',2^(1-nref), xi, omega));
%legend('PWDG L2','ultra-weak L2');

% LEGEND
leg1={};
leg2={};
for j=1:n_om
    leg11=strcat('\omega = ', num2str(omegas(j)) );
    leg1={leg1{1:j-1},leg11};
    leg22=strcat('PWDG, \omega = ', num2str(omegas(j)) );
    leg2={leg2{1:j-1},leg22};
end
legend(leg1{:}, leg2{:}, 'fontsize', 15)
%legend(leg1{:})
set(gca,'fontsize',18);
set(legend,'location','eastoutside')

% L2 - loglog
figure;
set(gca,'fontsize',18);
%loglog(ndir_l,err_uwvf,'--',ndir_l,err_pwdg,'-.');
%loglog(ndir_l,err_uwvf,'-','linewidth',2);
for j=1:n_om
    loglog(ndir_l,err_uwvf(j,:),'-', 'color',[0,0,1-(j-1)/(n_om-1)],'linewidth',2);
    hold on
end
grid off;
xlabel('\bf p/log(p)','fontsize',18);
ylabel('\bf L^2 error','fontsize',18);
%title(sprintf('%s Convergence of PWDG, h = %g, xi = %d, omega = %d','\bf',2^(1-nref), xi, omega));
%legend('PWDG L2','ultra-weak L2');
set(gca,'fontsize',18);

% LEGEND
leg1={};
leg2={};
for j=1:n_om
    leg11=strcat('\omega = ', num2str(omegas(j)) );
    leg1={leg1{1:j-1},leg11};
    leg22=strcat('PWDG, \omega = ', num2str(omegas(j)) );
    leg2={leg2{1:j-1},leg22};
end
legend(leg1{:}, leg2{:}, 'fontsize', 15)
%legend(leg1{:})

set(gca,'XLim',[min(ndir./log(ndir)),max(ndir./log(ndir))]);
legend off



%plot_Mesh(Mesh);   axis on

toc
