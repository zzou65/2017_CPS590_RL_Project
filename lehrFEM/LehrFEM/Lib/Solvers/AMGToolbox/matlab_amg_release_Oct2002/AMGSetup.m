%AMGSetup	AMG setup phase
%
%usage : [L,times] = AMGSetup( A [[,options] ,grid] )
%
%A       = matrix, matrix filename or matrix number in Names cell array
%options = options structure (defaults set by AMGDefaultOptions) or number
%          to indicate an option in the Opt struct array from a loaded 
%          Results file. Use a cell array to specify different options
%          for different levels: {Level1Options,Level2Options,...}. The
%          last option is used for the remaining levels.
%grid    = grid data file name, defaults to matrix file name - .rua + .grid,
%          if the matrix file name is given
%L       = returned data structure
%times   = structure with timing info


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These programs were prepared by the Regents of the University of
% California at Los Alamos National Laboratory (the University) under
% contract W-7405-ENG-36 with the U.S. Department of Energy (DOE) and
% under Contract KC-07-01-01 with the Department of Energy(DOE),
% Office of Science, Mathematical, Information, and Computational Sciences,
% Applied Mathematical Sciences Program.  All rights in these programs
% are reserved by the DOE and the University.

% Permission is granted to the public to copy and use this software
% without charge, provided that this Notice and the statements of 
% authorship are reproduced on all copies. Neither the U.S. government 
% nor the University makes any warranty, express or implied, or assumes 
% any liability or responsibility for the use of this software.

% AMG code written by Menno Verbeek, in collaboration with Jane Cullum
%  and Wayne Joubert.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function [L,times]=AMGSetup(A,options,grid)

%% Initialize

if nargin<3
  grid=[];
end
if length(A)==1
  disp(['Using name ',num2str(A),' from Names cell array'])
  A=evalin('base',['Names{',num2str(A),'}']);
end
if ischar(A)
  disp(['Reading matrix file ',A]);
  if isempty(grid)
    grid=strrep(A,'.rua','.grid');
  end
  A=readhb(A);
end

if nargin<2
  options=[];
end
if isnumeric(options) && length(options)==1
  disp(['Using option ',num2str(options),' from Opt struct array'])
  options=evalin('base',['Opt(',num2str(options),')']);
end
if iscell(options) && length(options)==1
  options=options{1};
end

if iscell(options)
  opt=options{1};
else
  opt=options;
end
DefaultOptions=AMGDefaultOptions;
opt=mergestruct(DefaultOptions,opt);
if ~opt.R.PT
  opt.R=mergestruct(DefaultOptions.P,opt.R);
end
L{1}.opt=opt;

inittimer;
starttimer('Total_Setup');

%% Do grid stuff

if isfield(opt,'mesh') && isfield(opt.mesh,'mesh')
  [grid,L{1}.meshx,L{1}.meshy]=mesh2coord(opt.mesh.mesh,opt.mesh.m,opt.mesh.n);
  if isfield(opt.mesh,'grid')
    L{1}.grid=opt.mesh.grid;
  else
    L{1}.grid=grid;
  end
elseif isfield(opt,'mesh') && isfield(opt.mesh,'grid')
  L{1}.grid=opt.mesh.grid;
  L{1}.meshx=[];
  L{1}.meshy=[];
elseif length(grid)>0
  D=dir(grid);
  if length(D)>1
    disp(['Warning, multiple grid files, using ',D(1).name]);
  end
  if length(D)>0
    L{1}.mesh=load(D(1).name);
    [t,n]=stripnr(D(1).name);
    t=t(1:end-1);
    [t,m]=stripnr(t);
    if n>0 && m>0
      disp(['Loading mesh size ',num2str(m),'x',num2str(n)]);
      [L{1}.grid,L{1}.meshx,L{1}.meshy]=mesh2coord(L{1}.mesh,m,n);
    else
      disp(['Could not get mesh size from file name ',D(1).name]);
    end
  end
else
  disp('No grid data')
end

%% Construct fine levels

l=0;
done=0;

while ~done
  l=l+1;

  if iscell(options) && length(options)>=l
    opt=mergestruct(opt,options{l});
  end
  L{l}.opt=opt;

  disp(['Setup of level ',num2str(l)]);
  if l==1
    L{l}.A=A;
  else
    L{l}.A=L{l-1}.Ac;
  end

  Symmetric = norm(A-A.',1)<1e-10*norm(A,1);

%% Coarse gird selection
  disp(' Coarse grid selection')
  starttimer('Coarse_grid_selection');
  [L{l}.nc,L{l}.perm,L{l}.Astrong,L{l}.type]=...
                                 AMGSelectCoarseGrid(L{l}.A,opt.CF);
  [t,f]=stoptimer('Coarse_grid_selection');
  disp(['   Time = ',num2str(t),' MFlops = ',num2str(f/1e6)]);

  if ~Symmetric && 0
    disp(' Coarse grid selection for transpose')
    starttimer('Coarse_grid_selection');
    [L{l}.nct,L{l}.permt,L{l}.Atstrong,L{l}.typet]=...
                                   AMGSelectCoarseGrid(L{l}.A',opt.CF);
    [t,f]=stoptimer('Coarse_grid_selection');
    disp(['   Time = ',num2str(t),' MFlops = ',num2str(f/1e6)]);
  else
    L{l}.nct=L{l}.nc;
    L{l}.permt=L{l}.perm;
    L{l}.Atstrong=L{l}.Astrong;
    L{l}.typet=L{l}.type;
  end  

  % Coarse-Fine points administration tools
  L{l}.c=L{l}.perm(1:L{l}.nc);
  L{l}.f=L{l}.perm(L{l}.nc+1:end);
  [dummy,L{l}.invperm]=sort(L{l}.perm);

  % Grid data
  if isfield(L{l},'grid')
    L{l+1}.grid=L{l}.grid(L{l}.c,:);
  end

%% Smoothers
  disp(' Pre-Smoother')
  starttimer('Pre_Smoother');
  [L{l}.Mpre,L{l}.Mpre_explicit]=GetSmoother(L{l},opt.pre);
  [t,f]=stoptimer('Pre_Smoother');
  disp(['   Time = ',num2str(t),' MFlops = ',num2str(f/1e6)]);
  disp(' Post-Smoother')
  starttimer('Post_Smoother');
  if opt.post.pre || comparestruct( opt.post, opt.pre )
    disp('   Mpost = Mpre');
    L{l}.Mpost=L{l}.Mpre;
    L{l}.Mpost_explicit=L{l}.Mpre_explicit;
  else
    [L{l}.Mpost,L{l}.Mpost_explicit]=GetSmoother(L{l},opt.post);
  end
  [t,f]=stoptimer('Post_Smoother');
  disp(['   Time = ',num2str(t),' MFlops = ',num2str(f/1e6)]);

%% Interpolation and restriction
  disp(' Interpolation')
  starttimer('Interpolation');
  L{l}.P=AMGMakeInterpolation( L{l}, opt.P, 0 );
  if exist('smooth','var')%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% added ",'var'"
    L{l}.smooth=smooth;
  end
  [t,f]=stoptimer('Interpolation');
  disp(['   Time = ',num2str(t),' MFlops = ',num2str(f/1e6)]);
  disp(' Restriction')
  starttimer('Restriction');
  if opt.R.PT || ( comparestruct( opt.P, opt.R ) && Symmetric )
    disp('   R = P^T');
    L{l}.R=L{l}.P.';
  else
    L{l}.R=AMGMakeInterpolation( L{l}, opt.R, 1 );
  end
  [t,f]=stoptimer('Restriction');
  disp(['   Time = ',num2str(t),' MFlops = ',num2str(f/1e6)]);

% jane: added nov 17, 1999 to use with Dff correction
  if opt.P.mDff == 0
% Not using Dff correction
    disp([' Not Using Dff correction: Level=',num2str(l)]);
    L{l}.imDff = 0;
%   L{l}.mDff = sparse(length(L{l}.perm),length(L{l}.perm));
  elseif ((opt.P.mDff == -1)&&(strcmp(opt.P.method,'R-S') || strcmp(opt.P.method,'Ruge-Stueben')))
% Using Dff-tilde correction: Only do if also doing R-S
    disp([' Using R-S and Dff-tilde Correction: Level=',num2str(l)]);
    L{l}.imDff = 1;
    [L{l}.mDff]=AMGMakeMDff(L{l},opt.P);
  elseif(opt.P.mDff == 1)
    disp([' Using Dff Correction from A: Level=',num2str(l)]);
    L{l}.imDff = 1;
    L{l}.mDff = spdiag(diag(L{l}.A(L{l}.f,L{l}.f)));
  else
    disp([' Unrecognizable opt.P.mDff option: Not Using Dff',num2str(l)]);
    L{l}.imDff = 0;
%   L{l}.mDff = sparse(length(L{l}.perm,L{l}.perm));
  end;
% jane: added to use Exact AFF-inverse in Dff Correction: if using Exact method
  if strcmp(opt.P.method,'exact')
    if (opt.P.exact == 1 && opt.P.mDff == 1)
% use Dff inverse correction from A
     L{l}.exact =1;
     L{l}.imDff = 1;
disp([' Exact but Using Dff Correction from A: Level=',num2str(l)]);
    elseif(opt.P.exact == -1 && opt.P.mDff == 1)
% use exact Aff-inverse correction
disp([' Exact and  Using Aff Correction from A: Level=',num2str(l)]);
     L{l}.exact = -1;
     L{l}.imDff = 1;
    else
% we should not be here: 
disp([' Exact: Not doing Dff Correction: Level=',num2str(l)]);
     L{l}.exact = 0;
     L{l}.imDff = 0;
    end;
  else
    L{l}.exact = 0;
  end;
% see below for AMGMakeMDff subroutine addition

%% Coarse grid matrix
  disp(' Coarse grid matrix')
  starttimer('Coarse_grid_matrix');
  L{l}.Ac=L{l}.R*L{l}.A*L{l}.P;
  [t,f]=stoptimer('Coarse_grid_matrix');
  disp(['   Time = ',num2str(t),' MFlops = ',num2str(f/1e6)]);

  if opt.levsmax==0
    done=( L{l}.nc <= opt.mincoarse );
  else
    done=( l == (opt.levsmax-1) );
  end
  if L{l}.nc/length(L{l}.A) > opt.max_coarse_frac
    % Don't do this last level
    disp(['Warning : coarse fraction too large, nc/n = ',...
          num2str(L{l}.nc/length(L{l}.A))]);
    L=L(1:l-1);
    l=l-1;
    % Grid data
    if isfield(L{l},grid)
      L{l+1}.grid=L{l}.grid(L{l}.c,:);
    end
    done=1;
  end
% jane: nov 17, 1999 added check if nc=0 to prevent bombing
  if L{l}.nc == 0
    % Don't do this last level
    disp(['Warning : nc = 0: level = ',...
          num2str(l)]);
    L=L(1:l-1);
    l=l-1;
    % Grid data
    if isfield(L{l},grid)
      L{l+1}.grid=L{l}.grid(L{l}.c,:);
    end
    done=1;
    error([' nc = 0: level= ',num2str(l)]);
% end of addition
  end
% end of while ~done, ie end of loop over fine levels
end

%% Do bottom level

l=l+1;
disp(['Setup of final level ',num2str(l)]);
if l==1
  error('Final level is first level => would be direct method\n');
  %L{l}.A=A;
else
  L{l}.A=L{l-1}.Ac;
end

disp(' LU-decomosition of A')
starttimer('LU');
[L{l}.L,L{l}.U]=lu(L{l}.A);
[t,f]=stoptimer('LU');
disp(['   Time = ',num2str(t),' MFlops = ',num2str(f/1e6)]);

[t,f]=stoptimer('Total_Setup');
disp(['Total Setup : Time = ',num2str(t),' MFlops = ',num2str(f/1e6)]);


times=getalltimers;

return


%% %%%%%%% Dff Corrections
% jane: added nov 17 1999
function[Dff]= AMGMakeMDff(L,opt)

if strcmp(opt.weak,'lump')
   f=L.f;
   Dff = spdiag(sum(L.A(f,:).') - sum(L.Astrong(f,:).'));
else
   Dff = spdiag(diag(L.A(f,f)));
end;
return


%% %%%%%%% Smoothers

function [M,explicit]=GetSmoother(Ll,opt)

  n=size(Ll.A,1);

  if isfield(opt,'M')
    disp('   Using given smoother (opt.M)');
    M=opt.M;
    explicit=opt.explicit;
  elseif strcmp( opt.type, 'Jacobi' )
    disp(['   Jacobi smoothing, omega=',num2str(opt.jac_omega)]);
    % Jacobi smoothing
% jane: Changed Nov 1 1999: 
%    M=JacobiSmoother(Ll.A);
% added nov 5 1999: opt.jac_omega=0 means compute upper bound on lambda(D^{-1}A)
%  opt.jac_omega ~0 means use the supplied value
    M=JacobiSmoother(Ll.A,opt.jac_omega);
    explicit=1;
  elseif strcmp( opt.type, 'Gauss-Seidel CF' )
    disp('   CF Gauss-Seidel smooting');
    % CF Gauss-Seidel smoothing    
    if isfield(opt,'GSperm') && length(opt.GSperm)==n
      disp('   Using permutation from opt.GSperm');
      GSperm=opt.GSperm;
      % Check whether we still have CF ordering
      if length(find(Ll.type(GSperm(1:Ll.nc    ))~=1))~=0 || ...
         length(find(Ll.type(GSperm(Ll.nc+1:end))~=2))~=0 
         disp('Warning : opt.GSperm is not in CF ordering');
      end
    else
      disp('   Warning : No input permutation in opt (opt.GSperm)');
      GSperm=Ll.perm;
    end
    M=GaussSeidelSmoother(Ll.A,GSperm);
    explicit=0;

  elseif strcmp( opt.type, 'Gauss-Seidel FC' )
    disp('   FC Gauss-Seidel smooting');
    % FC Gauss-Seidel smoothing    
    if isfield(opt,'GSperm') && length(opt.GSperm)==n
      disp('   Using permutation from opt.GSperm');
      GSperm=opt.GSperm;
      % Check whether we still have FC ordering
      if length(find(Ll.type(GSperm(1:end-Ll.nc-1))~=2))~=0 || ...
         length(find(Ll.type(GSperm(end-Ll.nc:end))~=1))~=0 
         disp('Warning : opt.GSperm is not in FC ordering');
      end
    else
      disp('   Warning : No input permutation in opt.GSperm');
      GSperm=Ll.perm;
    end
    M=GaussSeidelSmoother(Ll.A,GSperm);
    explicit=0;

  elseif strcmp( opt.type, 'Gauss-Seidel upper CF' )
    disp('   upper CF Gauss-Seidel smooting');
    % upper CF Gauss-Seidel smoothing    
    if isfield(opt,'GSperm') && length(opt.GSperm)==n
      disp('   Using permutation from opt.GSperm');
      GSperm=opt.GSperm;
      % Check whether we still have CF ordering
      if length(find(Ll.type(GSperm(1:Ll.nc    ))~=1))~=0 || ...
         length(find(Ll.type(GSperm(Ll.nc+1:end))~=2))~=0 
         disp('Warning : opt.GSperm is not in CF ordering');
      end
    else
      disp('   Warning : No input permutation in opt (opt.GSperm)');
      GSperm=Ll.perm;
    end
    M=GaussSeidelSmoother(Ll.A',GSperm)';
    explicit=0;

  elseif strcmp( opt.type, 'ainv' )
    disp(['Using ainv for smoother, tol=',num2str(opt.ainv_tol),...
          ' omega=',num2str(opt.ainv_omega),...
          ' shift=',num2str(opt.ainv_shift)]);
    [Z,D]=ainv(Ll.A+opt.ainv_shift*speye(size(Ll.A)),...
	       1,opt.ainv_tol);
    M=opt.ainv_omega*Z*inv(D)*Z.';
    explicit=1;

  else
    error(['   Unknown smoother type = ',opt.type]);
  end

return


%%%%%
% modified 11/01/99: Jane
% Needs explicit option
function M=JacobiSmoother(A,jomega )
  M=spdiag(1.0./diag(A));
  if jomega == 0
% Compute upper bound on Max(lambda(B)) where B=D^{-1}*A))
% Gershgorin Circles.
   rsB= (sum(abs(A'))'-abs(diag(A)))./abs(diag(A));
   clear omega;
   omega=max(1+rsB);
% Modify M: Omega is upper bound on Max(lambda(B)
   M=M/omega;
  else
% use jomega
   M=jomega*M;
  end;
return

%%%%%

function M=GaussSeidelSmoother(A,GSperm)
  [dummy,GSinvperm]=sort(GSperm);
  %M=inv(tril(A(GSperm,GSperm),0));
  M=tril(A(GSperm,GSperm),0);
% perhaps find smoother paramter based on
% power method or arnoldi
  iparm = 0;
  if iparm ~= 0
    rand('state',728591);
    v=rand(size(A,2),1);
    nr=6;
    nt=1;
    while nt <= nr
      v=M\(A*v);
      v=v/norm(v);
    end;
    nt=nt+1;
    omega= v'*(M\(A*v));
    M=omega*M;
  end;
  M=M(GSinvperm,GSinvperm);
return

%%%%%

function [name,nr]=stripnr(name)
  name_=[name,' '];
  numbers=(name_>='0')&(name_<='9');
  start=find(numbers(1:end-1)==0 & numbers(2:end)==1)+1;
  stop =find(numbers(1:end-1)==1 & numbers(2:end)==0);
  if length(start)>0
    nr=str2num(name(start(end):stop(end)));
    name=name(1:start(end)-1);
    name_=find(name~='_');
    name=name(1:name_(end));
  else
    nr=-1;
  end
return


