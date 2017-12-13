%AMGDefaultOptions	Return default options structure
%
%usage : opt=AMGDefaultOptions


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% These programs were prepared by the Regents of the University of
%% California at Los Alamos National Laboratory (the University) under
%% contract W-7405-ENG-36 with the U.S. Department of Energy (DOE) and
%% under Contract KC-07-01-01 with the Department of Energy(DOE),
%% Office of Science, Mathematical, Information, and Computational Sciences,
%% Applied Mathematical Sciences Program.  All rights in these programs
%% are reserved by the DOE and the University.

%% Permission is granted to the public to copy and use this software
%% without charge, provided that this Notice and the statements of 
%% authorship are reproduced on all copies. Neither the U.S. government 
%% nor the University makes any warranty, express or implied, or assumes 
%% any liability or responsibility for the use of this software.

%% AMG code written by Menno Verbeek, in collaboration with Jane Cullum
%%  and Wayne Joubert.
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function opt=AMGDefaultOptions

% Levels
opt.levsmax=0;
opt.mincoarse=100;
opt.max_coarse_frac=0.75;

% Pre-Smoother
opt.pre.its=2;
opt.pre.type='Gauss-Seidel CF';  % 'Jacobi', 'Gauss-Seidel CF', 
                                 % 'Gauss-Seidel upper CF',
                                 % 'Gauss-Seidel FC', 'ainv'
opt.pre.GSperm=[];
% jane: opt.pre.jac_omega=0: compute and use upper bound on
%   maximum eigenvalue of D^{-1}*A: opt.pre.jac_omega ~= 0 means
%   multiply D^{-1}*A by that value.
% changed omega to 0 on 3/27/2001. 1 option does not
%  work well.

opt.pre.jac_omega=0;
opt.pre.ainv_tol=0.1;
opt.pre.ainv_omega=1;
opt.pre.ainv_shift=0.05;
%opt.pre.M=[...];

% Post-Smoother
opt.post.its=2;
opt.post.pre=1;

% Coarse grid selection
opt.CF.theta=0.25;
opt.CF.use_abs=0;
opt.CF.use_diag=0;
opt.CF.direction='row';           % 'row', 'row and column'
% changed on Nov 19 1999 from 0 to 1
opt.CF.make_unconnected_F=1;
opt.CF.check_F_F=1;
%opt.CF.perm=[...];
%opt.CF.nc=..;

% Interpolation
opt.P.method    = 'Ruge-Stueben'; % 'Ruge-Stueben'/'R-S', 'exact', 'ainv',


                                  % 'ainv strong', 'gsai', 'local','locallocal'
opt.P.pattern   = 'Astrongfc';    % 'dense', 'Afc', 'Astrongfc'
opt.P.matrix    = 'A';            % 'A', 'MA', 'Se', 'Sr'
opt.P.objective = '-Aff\Afc';     % '0', '-Aff\Afc', '-A:f\A:c', 'Aff/Acf',
                                  % 'eigenvecs', 'Left-singular-vecs',
                                  % 'Right-singular-vecs'
opt.P.droptol   = 0;
opt.P.ainv_tol  = 0.1;
opt.P.weak      = 'lump';         % 'lump', 'scale', 'trash'
opt.P.method64  = 'original';     % 'original', 'original abs', 'original abs abs'
                                  % 'new1', 'new2'
opt.P.nofill64  = 1;              % 0, 1
opt.P.eig_use_distance = 0;

%jane: added Nov 17 1999
% Allows adding correction to basic R-S AMG Fine x correction
opt.P.mDff = 0;        % 0 = do not use; 1 = use; when using ainv set to 0
% =0 if method ~= 'exact':  =1 if method ='exact' and opt.P.imDff=1; =0 otherwise.
opt.P.exact = 0;

% Restriction
opt.R.PT = 1;

% Mesh
%opt.mesh.mesh=[...];
%opt.mesh.m=..;
%opt.mesh.n=..;
%opt.mesh.grid=[...];

