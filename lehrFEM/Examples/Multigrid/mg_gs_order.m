function [] = mg_gs_order(smoother)
% effects of ordering of smoothing iterations
%
%   This code plots the energy-norm error vs the iteration number for
%   multigrid V-cycles with different smoothing orders,
%     CF-CF : smooth coarse nodes first in presmoothing and postsmoothing
%     FC-FC : smooth fine nodes first in presmoothing and postsmoothing
%     CF-FC : smooth coarse nodes first in presmoothing and fine nodes
%             first in postsmoothing
%     FC-CF : smooth fine nodes first in presmoothing and coarse nodes
%             first in postsmoothing
%     random : use random permutation
%
%   In the sample problem, the nonsymmetric orderings CF-CF and FC-FC
%   work best for the Gauss-Seidel smoother but the symmetric orderings
%   CF-FC and FC-CF work best for the symmetric Gauss-Seidel smoother (use
%   argument @sgs_smooth).
%
%   If symmetric Gauss-Seidel is interpereted as two Gauss-Seidel steps,
%   then the methods that work best are
%     CF-CF, FC-FC, CFFC-FCCF, FCCF-CFFC
%   and those that converge somewhat slower are
%     CF-FC, FC-CF, CFFC-CFFC, FCCF-FCCF (and random).
%   It seems that the 'skew-symmetric' methods are more efficient than the
%   symmetric methods.

%   Copyright 2007-2007 Claude Gittelson
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  ref = [0 7];                                   % refinements of mesh
  f_handle = @(x,varargin)-4*ones(size(x,1),1);  % Right hand-side source term
  gd_handle = @(x,varargin)x(:,1).^2+x(:,2).^2;  % Dirichlet boundary data
%   f_handle = @f_LShap;                           % Right hand-side source term
%   gd_handle = @g_D_LShap;                        % Dirichlet boundary data

  if(nargin==0)
    smoother = @gs_smooth;                         % multigrid smoother
%   smoother = @sgs_smooth;                        % multigrid smoother
  end
  cyc = 1;                                       % V or W multigrid cycles
  m = [1 1];                                     % number of smoothing steps
  
  tol = 1e-6;                                    % Stopping criterion
  maxit = 100;                                   % Maximum number of iterations
  
  % generate multigrid data structure
  
  CMesh = load_Mesh('Coord_Sqr.dat','Elem_Sqr.dat');
%   CMesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
  mg_data = mg_mesh('mesh',CMesh,'ref',ref);
  mg_data = mg_stima(mg_data,'f',f_handle,'gd',gd_handle);
  mg_data = mg_error(mg_data,'iter',false,'exact',true,'rel',false);
  LVL = length(mg_data);
  
  % initialize data
  
  x0 = zeros(size(mg_data{LVL}.b));
  err0 = mg_data{LVL}.error.energy_exact(x0);
  
  % run multigrid solver with coarse nodes smoothed first
  
  mg_data = mg_smooth(mg_data,'m',m,'cyc',cyc,'smoother',smoother,...
    'per','CF','sym_per',false);
  [x,conv] = mg(x0,mg_data,tol,maxit);
  numit_CF = conv.iter;
  err_CF = [err0,conv.error.energy_exact];
  
  % run multigrid solver with coarse nodes smoothed last
  
  mg_data = mg_smooth(mg_data,'m',m,'cyc',cyc,'smoother',smoother,...
    'per','FC','sym_per',false);
  [x,conv] = mg(x0,mg_data,tol,maxit);
  numit_FC = conv.iter;
  err_FC = [err0,conv.error.energy_exact];
  
  % run multigrid solver with coarse nodes smoothed first in presmoothing
  % and last in postsmoothing
  
  mg_data = mg_smooth(mg_data,'m',m,'cyc',cyc,'smoother',smoother,...
    'per','CF','postper','FC','sym_per',false);
  [x,conv] = mg(x0,mg_data,tol,maxit);
  numit_CF_FC = conv.iter;
  err_CF_FC = [err0,conv.error.energy_exact];
  
  % run multigrid solver with coarse nodes smoothed last in presmoothing and
  % first in postsmoothing
  
  mg_data = mg_smooth(mg_data,'m',m,'cyc',cyc,'smoother',smoother,...
    'per','FC','postper','CF','sym_per',false);
  [x,conv] = mg(x0,mg_data,tol,maxit);
  numit_FC_CF = conv.iter;
  err_FC_CF = [err0,conv.error.energy_exact];
  
  % run multigrid solver with random order in smoother
  
  mg_data = mg_smooth(mg_data,'m',m,'cyc',cyc,'smoother',smoother,...
    'per','rand','sym_per',false);
  [x,conv] = mg(x0,mg_data,tol,maxit);
  numit_rand = conv.iter;
  err_rand = [err0,conv.error.energy_exact];
  
  % plot error vs. number of v-cycle iterations
  
  figure;
  semilogy(0:numit_CF,err_CF,'-o',...
           0:numit_FC,err_FC,'-s',...
           0:numit_CF_FC,err_CF_FC,'-<',...
           0:numit_FC_CF,err_FC_CF,'->',...
           0:numit_rand,err_rand,'-+');
         
  legend('CF-CF','FC-FC','CF-FC','FC-CF','random','Location','NorthEast');
  grid on;
  set(gca,'YLim',[tol,err0]);
  xlabel('\bf V-cycles');
  ylabel('\bf error in energy norm');
  title('\bf Effects of Ordering of the Smoothing Iterations');

return