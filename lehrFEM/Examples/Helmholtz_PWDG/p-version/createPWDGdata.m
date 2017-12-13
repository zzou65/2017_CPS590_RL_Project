% 01.06.2010 - Andrea - Barcelona
% solve Helmholtz on a square with different parameters
% save errors and parameters in PWDG_error_data.mat
% for various plots and for Jay


clear all; close all;format compact;

% parameters
ndir=[3:2:41];
omegas=[0.125 0.25 0.5 1 2 4 8 16 32 64];
nref=[1 2 3];
xis=[1 2/3 3/2];
MyFileName = 'PWDG_error_data';
UseFile = 1;


%TEMPORAY PARAMETERS....................
%   ndir=[3:2:15];
%   omegas=[ 8 10];
%   nref=[1];
%   xis=[1];

%load file if it exists and UseFile is on, PWDG will contain all the data
if UseFile 
    try
        load (MyFileName);
    catch
        PWDG = {};
    end
else
    PWDG = {};
end
ind = length(PWDG); %index in PWDG{}


%% Define discretization parameters             

% Define flux parameters
flcoef=10;  %coefficient to improve the convergence of the method with PWDG fluxes
flux_params_pwdg = {'a',@(omega,h,p,varargin) flcoef*p/(h*omega*log(p)),   'b',@(omega,h,p,varargin) h*omega*log(p)/(p*flcoef),'d',@(omega,h,p,varargin) h*omega*log(p)/(p*flcoef)};
flux_params_uwvf = {'a',0.5,'b',0.5,'d',0.5};

% Define quadrature rules (then I can decide where to use qr2 and qr3)
qr1 = gauleg(0,1,12);                   % 1D for edges
qr2 = Duffy(TProd(gauleg(0,1,12)));     % for triangular elements
qr3 = Duffy(TProd(gauleg(0,1,50)));     % for triangular elements  high order to resolve the singularities
% qr2 = TProd(gauleg(0,1,12)); % for quadrilateral elements

% square domain
R=1;
offsetR=0;                                                                                    % offset to move the mesh in order to have the singularity on a edge, not on a corner
Mesh.Coordinates = [0 -R/2 + offsetR; R -R/2 + offsetR; R R/2 + offsetR; 0 R/2 + offsetR];    % right square
Mesh.Elements = [1 2 3; 1 4 3]; % triangular elements
Mesh = add_Edges(Mesh);         
Loc = get_BdEdges(Mesh);
Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
Mesh.BdFlags(Loc) = -1;

%%%%%%% now the 4 loops in the parameters: [ref xi omega p]
%% LOOP IN MESHES
        
for ref=nref
    Mesh = refine_REG(Mesh);
    %plot_Mesh(Mesh)
    meshwidth = 2* R / sqrt(length(Mesh.Elements));    % only for regular triangular mesh on square RxR!
    
    Mesh = orient_Elems(Mesh);
    Mesh = add_Edge2Elem(Mesh);
    Mesh = add_DGData(Mesh);
    
    %% LOOP XI - OMEGA    
    for xi = xis
        for omega =omegas
            lambda=omega; %lambda is for rhs
            omegah= omega*meshwidth;
            % exact solutions
            vnorm=@(x) sqrt(x(:,1).^2+x(:,2).^2) ;          %euclid norm 2D for arrays N x 2 (--> result Nx1)
            vangle = @(x) angle(x(:,1)+1i*x(:,2));
            %circular wave for  xi integer,
            %Helmholtz corner singularity in the origin for real xi non integer
            theta_0=0;
            u_ex = @(x,varargin) (besselj(xi, lambda*vnorm(x) )) .*cos(xi*(vangle(x)-theta_0));
            
            Dx_u_ex= @(x,varargin) cos(vangle(x)) .* cos(xi*(vangle(x)-theta_0))*lambda/2 .* (besselj(xi-1, lambda*vnorm(x)) - besselj(xi+1, lambda*vnorm(x) ) ) + ...
                xi* sin(vangle(x)) .* sin(xi*(vangle(x)-theta_0)) ./ vnorm(x) .* besselj(xi,  lambda*vnorm(x));
            Dy_u_ex= @(x,varargin) sin(vangle(x)) .* cos(xi*(vangle(x)-theta_0))*lambda/2 .* (besselj(xi-1, lambda*vnorm(x)) - besselj(xi+1, lambda*vnorm(x) ) ) - ...
                xi* cos(vangle(x)) .* sin(xi*(vangle(x)-theta_0)) ./ vnorm(x) .* besselj(xi,  lambda*vnorm(x));
            
            % impedence and Dirichlet condition if not already defined
            gI = @(x,n,varargin) [Dx_u_ex(x), Dy_u_ex(x)] * n'  + 1i*omega* u_ex(x);
            %gD=u_ex;
            f= @(x,varargin) 0;             %source term
            grad_u_ex= @(x,varargin) [Dx_u_ex(x,varargin), Dy_u_ex(x,varargin) ];
            
              
        %% LOOP IN p
            for p = ndir
                tic
            
            %now check if configuration [R xi omega ref p] already exists,
            %if so, it skips the computation            % from Roman
                skip = false;
            for ii = 1:length(PWDG)
                data = PWDG{ii};
                if ((data.R == R) && (data.xi == xi) && (data.omega == omega) && (data.ref == ref) && (data.p == p))
                    skip = true;
                    break;
                end
            end
            if (skip);        disp('Already available');        continue;    end
            
            % keep ind+1 after the check on the file... sure?
            ind = ind + 1
            
            %# degrees of freedom =  size lin. system
            ndof = p * length(Mesh.Elements);
            disp('===================================');
            disp([' #refinements = ' num2str(ref)]);
            disp([' xi = ' num2str(xi)]);
            disp([' omega = ' num2str(omega)]);
            disp([' h = ' num2str(meshwidth)]);
            disp([' p = ' num2str(p)]);
            disp([' ndof = ' num2str(ndof)]);
            disp('-----------------------------------');
            %Mesh
            
  %now build 3 matrices, 3 rhs and solve 3 systems
                                  
  % Define plane wave directions
  Mesh = set_Data_PWDG(Mesh,'nDofs',p,'NewDir','replace','Omega',omega);
  % Set flux parameters for plane wave DG
  Mesh = set_Data_PWDG(Mesh,'Omega',omega,flux_params_pwdg{:});
  % Assemble stiffness matrix for plane waves
  [I_Inn,J_Inn,A_Inn] = assemMat_Inn_PDG2(Mesh,@STIMA_Helm_Inn_PWDG,omega);
      [I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2(Mesh,[],@STIMA_Helm_Imp_Bnd_PWDG,omega);
  A_pwdg = sparse([I_Inn;I_Bnd],[J_Inn;J_Bnd],[A_Inn;A_Bnd]);
  %condition number:
    try  cond_pwdg = condest(A_pwdg);
    catch; cond_pwdg = nan; end;

  % Assemble load vector for plane waves
    b_g = assemLoad_Bnd_PDG2(Mesh,[],@LOAD_Imp_Bnd_PWDG,qr1,omega,gI);
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
  A_uwvf = sparse([I_Inn;I_Bnd],[J_Inn;J_Bnd],[A_Inn;A_Bnd]);
  %condition number:
    try  cond_uwvf = condest(A_uwvf);
    catch;cond_uwvf = nan;  end;

 
  % Assemble load vector for ultraweak variational formulation
    b_g = assemLoad_Bnd_PDG2(Mesh,[],@LOAD_Imp_Bnd_PWDG,qr1,omega,gI);
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
  %condition number:
    try  cond_proj = condest(M);
    catch;cond_proj = nan;  end;
  
  % Calculate L2 projection of exact solution onto discretization space
  b_proj = assemLoad_Vol_PDG2(Mesh,@LOAD_Vol_PWDG,qr2,omega,u_ex);
  u_proj = M\b_proj;
 
  [err_pwdg_L2, err_pwdg_h1, err_pwdg_dg] = MultiNormErr_PWDG(Mesh,u_pwdg,qr2,omega,u_ex, grad_u_ex    ,1);           %quadratures! ...
  [err_uwvf_L2, err_uwvf_h1, err_uwvf_dg] = MultiNormErr_PWDG(Mesh,u_uwvf,qr2,omega,u_ex, grad_u_ex    ,1);
  [err_proj_L2, err_proj_h1, err_proj_dg] = MultiNormErr_PWDG(Mesh,u_proj,qr2,omega,u_ex, grad_u_ex    ,1);
  
  err_pwdg_L2 = sqrt(err_pwdg_L2);
  err_pwdg_h1 = sqrt(err_pwdg_h1);
  err_pwdg_dg = sqrt(err_pwdg_dg);
  err_uwvf_L2 = sqrt(err_uwvf_L2);
  err_uwvf_h1 = sqrt(err_uwvf_h1);
  err_uwvf_dg = sqrt(err_uwvf_dg);
  err_proj_L2 = sqrt(err_proj_L2);
  err_proj_h1 = sqrt(err_proj_h1);
  err_proj_dg = sqrt(err_proj_dg);
  
  PWDG{ind}.p = p;
  PWDG{ind}.omega = omega;
  PWDG{ind}.R = R;
  PWDG{ind}.ref = ref;
  PWDG{ind}.xi = xi;
  PWDG{ind}.omegah = omegah;
  PWDG{ind}.meshwidth = meshwidth;
  PWDG{ind}.ndof = ndof;
  
  PWDG{ind}.err_pwdg_L2 = err_pwdg_L2;
  PWDG{ind}.err_pwdg_h1 = err_pwdg_h1;
  PWDG{ind}.err_pwdg_dg = err_pwdg_dg;
  PWDG{ind}.err_uwvf_L2 = err_uwvf_L2;
  PWDG{ind}.err_uwvf_h1 = err_uwvf_h1;
  PWDG{ind}.err_uwvf_dg = err_uwvf_dg;
  PWDG{ind}.err_proj_L2 = err_proj_L2;
  PWDG{ind}.err_proj_h1 = err_proj_h1;
  PWDG{ind}.err_proj_dg = err_proj_dg;
  
  PWDG{ind}.cond_pwdg = cond_pwdg;
  PWDG{ind}.cond_uwvf = cond_uwvf;
  PWDG{ind}.cond_proj = cond_proj;
  
  if UseFile
      save (MyFileName, 'PWDG');
  end
                
%this toc contains total time for check into the file, assembly, solution,
%condition estimate, saving into the file...
                toc
           end
        end
    end
end

if UseFile
    save (MyFileName, 'PWDG');
    clear all;
    load PWDG_error_data
end


