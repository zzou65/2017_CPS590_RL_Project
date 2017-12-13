function [] = main4()
% convergence rates for discontinuous coefficients on line
%
%   Compares multigrid and algebraic multigrid for discontinuous
%   conductivity.  The conductivity is epsilon on the left half of the
%   square domain and one on the right half.  The boundary is resolved by
%   the coarsest mesh; both geomatric and algebraic multigrid perform well
%   for all epsilon.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
 
  REF = [2,6];                                   % refinements of initial mesh
  F_HANDLE = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  GD_HANDLE = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data
  M = [1,1];                                     % Smoothing steps
  
  MAXIT = 50;        % number of iterations in power method
  
  % Initialize multigrid data structure
  
  CMesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
  mg_data = mg_mesh('mesh',CMesh,'ref',REF);
  mg_data = mg_smooth(mg_data,'m',M);
  mg_data = mg_error(mg_data,'iter',false,'exact',false);
  
  Mesh = mg_data{1}.mesh;
  
  % Define epsilon
  
  epsilon = 10.^(0:-0.5:-4.5);
  c_mg = nan(size(epsilon));
  c_amg = nan(size(epsilon));
  
  for j=1:length(epsilon)

    % Define conductivity

    C_HANDLE = @(x,varargin) 1*(x(:,1)<0.5)...
                          + 0.5*(1+epsilon(j))*(x(:,1)==0.5)...
                          + epsilon(j)*(x(:,1)>0.5);

    % Update multigrid data structure

    mg_data = mg_stima(mg_data,'f',F_HANDLE,'gd',GD_HANDLE,'c',C_HANDLE);

    A = mg_data{end}.A;
    b = mg_data{end}.b;

    % Generate algebraic multigrid options

    AMGopt = AMGDefaultOptions;
    AMGopt.pre.its = M(1);
    AMGopt.post.its = M(end);
    AMGopt.mincoarse = min(AMGopt.mincoarse,mg_data{2}.n.free-1);

    % calculate convergence rate for multigrid using the power method;
    % the error propagation operator is E = I - BA
    
    Eerr = ones(size(b));
    u0 = zeros(size(b));
    for i=1:MAXIT
      err = Eerr;
      c = mg(u0,A*err,mg_data,0,1); % c = BAerr
      Eerr = err - c;
    end
    c_mg(j) = (err'*Eerr)/(err'*err);

    % calculate convergence rate for algebraic multigrid using the power method

    L = AMGSetup(A,AMGopt);

    Eerr = ones(size(b));
    for i=1:MAXIT
      err = Eerr;
      c = AMGVcycle(L,A*err);
      Eerr = err - c;
    end
    c_amg(j) = (err'*Eerr)/(err'*err);

  end % loop over epsilons
  
  % plot results
  
  figure;
  semilogx(epsilon,c_mg,'-o',epsilon,c_amg,'-+');
  xlabel('\bf \epsilon');
  ylabel('\bf maximal eigenvalue of error propagation matrix');
  title('\bf Convergence of Multigrid Methds for Discontinuous Coefficients on Line');
  legend('multigrid','algebraic multigrid','Location','NorthEast','Orientation','vertical');
  grid on;
  set(gca,'XDir','reverse','XLim',[min(epsilon),max(epsilon)],'YLim',[0,0.25]);
  
  % plot conductivity
  
  x = linspace(0,1);
  y = linspace(0,1);
  X = x(ones(size(y)),:);
  Y = y(ones(size(x)),:)';
  Z = zeros(size(X));
  Z(:) = C_HANDLE([X(:),Y(:)]);
  
  figure;
  contourf(x,y,Z);
  colorbar;
  plot_Mesh(Mesh,'f');
  
return  