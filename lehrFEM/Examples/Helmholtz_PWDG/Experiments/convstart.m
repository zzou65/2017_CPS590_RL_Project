function [fig1,fig2,homega,err_discr] = convstart(refs,omega_p,ndir,lambda_,d,flux_params,mtd_name,homega,err_discr)
%CONVSTART start of h-convergence
%   
%   CONVSTART(REFS,OMEGA_P,NDIR,LAMBDA_,D,PARAMS_INNER,PARAMS_BND,MTD_NAME)
%   plots the convergence of the plane wave discontinuous Galerkin method.
%
%   REFS is a vector of length 2 specifying the number of refinements on
%   the coarsest and finest mesh.
%
%   The wave number in the Helmholtz equation is 2.^OMEGA_P.  OMEGA_P
%   should be a vector of whole numbers.
%
%   NDIR is the number of local plane wave basis functions.
%
%   LAMBDA_ is a function handle with argument omega that determines the
%   frequency of the exact solution of the sample problem.
%
%   D is the direction of propagation of the exact solution.
%
%   FLUX_PARAMS is a cell array containing the flux parameter names and
%   values.  The default values are equivalent to
%       FLUX_PARAMS = {'a',0.5,'b',0.5,'c',[0 0],'d',0.5} .
%
%   See also run_convstart.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize discretization parameters
  nref = refs(2)-refs(1)+1;           % number of mesh refinements
  h = 2.^(1-refs(1):-1:1-refs(2));    % mesh width
  
  phi = 2*pi*(0:1/ndir:1-1/ndir)';    
  dir = [cos(phi),sin(phi)];          % directions of plane wave basis functions
  
  % Initialize problem parameters
  omega = 2.^omega_p;
  lambda = lambda_(omega);
  psi = 0;
  d = d/norm(d);
  nomega = numel(omega);
  
  % Define problem parameters
  u_ex = @(x,k,varargin) sin(lambda(k)*x*d'+psi);
  du_ex = @(x,k,varargin) lambda(k)*cos(lambda(k)*x*d([1 1],:)'+psi).*d(ones(size(x,1),1),:);
  g = @(x,n,k,varargin) lambda(k)*(d*n')*cos(lambda(k)*x*d'+psi) + i*omega(k)*sin(lambda(k)*x*d'+psi);
  f = @(x,k,varargin) (lambda(k)^2-omega(k)^2)*sin(lambda(k)*x*d'+psi);
  
  % Initialize quadrature rules
  qr1 = gauleg(0,1,20);
  qr2 = Duffy(TProd(gauleg(0,1,20)));
  
  % Initialize errors
  if(nargin<10)
    err_discr = nan(nomega,nref);
    homega_discr = nan(nomega,nref);
    do_discr = true;
  else
    do_discr = false;
  end
  err = nan(nomega,nref);
  
  % Initialize mesh
  Mesh.Coordinates = [-1 -1; 1 -1; 1 1; -1 1];
  Mesh.Elements = [1 2 3; 1 4 3];
  Mesh = add_Edges(Mesh);         
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
  Mesh.BdFlags(Loc) = -1;
  for j=1:refs(1)-1
    Mesh = refine_REG(Mesh);
  end
  
  % Loop over meshes
  for j=1:nref

    % Refine mesh and add fields
    Mesh = refine_REG(Mesh);
    Mesh = orient_Elems(Mesh);
    Mesh = add_Edge2Elem(Mesh);
    Mesh = add_DGData(Mesh);
    Mesh = set_Data_PWDG(Mesh,'Dir',dir);
    
    % Loop over values of omega
    for k=1:nomega
      
      % set flux parameters
      Mesh = set_Data_PWDG(Mesh,'Omega',omega(k),flux_params{:});
      
      %% Do discretization error if needed
      
      if(do_discr)
      
        % Assemble energy-norm inner product
        [I_h1_diag,J_h1_diag,H_diag] = assemMat_Vol_PDG2(Mesh,@SCAPRO_Helm_diag_PWDG,omega(k));
        [I_h1_inn,J_h1_inn,H_inn] = assemMat_Inn_PDG2(Mesh,@SCAPRO_Helm_Inn_PWDG,omega(k));
        [I_h1_bnd,J_h1_bnd,H_bnd] = assemMat_Bnd_PDG2(Mesh,[],@SCAPRO_Helm_Bnd_PWDG,omega(k));
        I_h1 = [I_h1_diag;I_h1_inn;I_h1_bnd];
        J_h1 = [J_h1_diag;J_h1_inn;J_h1_bnd];
        H_h1 = [H_diag;H_inn;H_bnd];
        B = sparse(I_h1,J_h1,H_h1);

        % Calculate projection of exact solution onto discretization space
        b_proj = assemLoad_Vol_PDG2(Mesh,@LOAD_EnergyProj_Vol_PWDG,qr2,omega(k),u_ex,du_ex,k);
        u_proj = B\b_proj;

        % Calculate error
        err_discr(k,j) = EnergyErr_PWDG(Mesh,u_proj,qr2,omega(k),u_ex,du_ex,0,k)/omega(k);
        homega_discr(k,j) = h(j)*omega(k);
        
      end
      
      %% Do full error
      
      % Assemble stiffness matrix
      [I_Inn,J_Inn,A_Inn] = assemMat_Inn_PDG2(Mesh,@STIMA_Helm_Inn_PWDG,omega(k));
      [I_Bnd,J_Bnd,A_Bnd] = assemMat_Bnd_PDG2(Mesh,[],@STIMA_Helm_Imp_Bnd_PWDG,omega(k));
      A = sparse([I_Inn;I_Bnd],[J_Inn;J_Bnd],[A_Inn;A_Bnd]);

      % Assemble load vector
      b_g = assemLoad_Bnd_PDG2(Mesh,[],@LOAD_Imp_Bnd_PWDG,qr1,omega(k),g,k);
      b_f = assemLoad_Vol_PDG2(Mesh,@LOAD_Vol_PWDG,qr2,omega(k),f,k);
      b = b_f + b_g;

      % Solve equation
      u = A\b;

      % Calculate error
      err(k,j) = EnergyErr_PWDG(Mesh,u,qr2,omega(k),u_ex,du_ex,0,k);
      
    end
  
  end
  
  % Extract discretization error
  if(do_discr)
    [homega,ind_discr] = unique(homega_discr(:));
    err_discr = err_discr(ind_discr);
  end
  
  % Plot errors with projection error
  fig1 = figure;
  lgd = cell(1,nomega);
  plot(homega,err_discr,'k--');
  hold on;
  for k=1:nomega
    plot(h*omega(k),err(k,:)/omega(k));
    hold all;
    lgd{k} = sprintf('%s = %g, %s = %g','\omega',omega(k),'\lambda',lambda(k));
  end
  hold off;
  set(gca,'XScale','log','XDir','reverse','YScale','log');
  grid on;
  legend('projection error',lgd{:},'Location','SouthWest');
  xlabel('\bf h\omega');
  ylabel('\bf ( \omega^{-2}| error |_{1,h}^2 + || error ||_0^2 )^{1/2}');
  title(sprintf('%s Convergence of %s for u(x) = sin(%sx%sd)','\bf',mtd_name,'\lambda','\cdot'));
  
  % Plot errors
  fig2 = figure;
  lgd = cell(1,nomega);
  for k=1:nomega
    plot(h*omega(k),err(k,:));
    hold all;
    lgd{k} = sprintf('%s = %g, %s = %g','\omega',omega(k),'\lambda',lambda(k));
  end
  hold off;
  set(gca,'XScale','log','XDir','reverse','YScale','log');
  grid on;
  legend(lgd{:},'Location','SouthWest');
  xlabel('\bf h\omega');
  ylabel('\bf ( | error |_{1,h}^2 + \omega^2|| error ||_0^2 )^{1/2}');
  title(sprintf('%s Convergence of %s for u(x) = sin(%sx%sd)','\bf',mtd_name,'\lambda','\cdot'));
  
return  