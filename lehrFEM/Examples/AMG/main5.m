function [] = main5()
% convergence rates for discontinuous coefficients on square
%
%   Compares multigrid and algebraic multigrid for discontinuous
%   conductivity.  The conductivity is epsilon inside a square in the
%   middle of the (also square) domain and one outside.  The inner square
%   is resolved by the coarsest mesh; both geometric and algebraic
%   multigrid perform well.

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
    
    c_handle = @(x,varargin) cond(x,epsilon(j));

    % Update multigrid data structure

    mg_data = mg_stima(mg_data,'f',F_HANDLE,'gd',GD_HANDLE,'c',c_handle);

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
  title('\bf Convergence of Multigrid Methds for Discontinuous Coefficients on Square');
  legend('multigrid','algebraic multigrid','Location','NorthEast','Orientation','vertical');
  grid on;
  set(gca,'XDir','reverse','XLim',[min(epsilon),max(epsilon)]);
  
  % plot conductivity
  
  x = linspace(0,1);
  y = linspace(0,1);
  X = x(ones(size(y)),:);
  Y = y(ones(size(x)),:)';
  Z = zeros(size(X));
  Z(:) = c_handle([X(:),Y(:)]);
  
  figure;
  contourf(x,y,Z);
  colorbar;
  plot_Mesh(Mesh,'f');
  
return


% the dependence of the conductivity on epsilon is given by the function

function  c = cond(x,epsilon)

  ind = x(:,1)>0.25 & x(:,1)<0.75 & x(:,2)>0.25 & x(:,2)<0.75;
  c = ones(size(x,1),1);
  c(ind) = epsilon;

return
  