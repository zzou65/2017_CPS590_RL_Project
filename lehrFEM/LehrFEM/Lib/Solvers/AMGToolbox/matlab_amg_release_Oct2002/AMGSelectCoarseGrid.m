%AMGSelectCoarseGrid	Selects the course grid points (Ruge-Stueben)
%
%usage : [nc,perm,Astrong,type]=AMGSelectCoarseGrid(A,opt)
%
%A   = the matrix
%opt = the CF part of the options (options.CF)
%nc  = the number of coarse points
%perm= the reordering vector, A(perm,perm)=[Acc,Acf]
%                                          [Afc,Aff]
%Astrong = the strong part of the matrix A, no diagonal
%type = vector with a 1 for coarse points, 2 for fine 
%       points (in original order) 


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


function [nc,perm,Astrong,type]=AMGSelectCoarseGrid(A,opt)

n=length(A);
U=3;
F=2;
C=1;


%%%% make Astrong

if ~opt.use_diag
  disp('   Not using diagonal of A for Astrong')
  Astrong=(A-spdiag(diag(A)));
else
  disp('   Using diagonal of A for Astrong')
  Astrong=A;
end

disp(['   Theta = ',num2str(opt.theta)])

if strcmp( opt.direction, 'row' )
  disp('   Using row criteria for Astrong');
  Keep=ThetaColumnCheck(Astrong.',opt).';
elseif strcmp( opt.direction, 'row and column' )
  disp('   Using row and column criteria for Astrong');
  Keep=ThetaColumnCheck(Astrong.',opt).';
  Keep=Keep | ThetaColumnCheck(Astrong,opt);
else
  error(['   Unknown direction for Astrong : ',opt.direction]);
end
Astrong=Keep.*A;
Astrong=(Astrong-spdiag(diag(Astrong)));


%%%%%%% Preform A2 and A3 in par. 6.2.6 from Runge-Stueben article

disp('   Do linked list version of first phase');

%%%% Set up connectivity numbers

lambda=full(sum(spones(Astrong)));

%%%% Put all points in the U set

type=U*ones(n,1);

%%%% A2 : Find the initial Coarse and Fine points

AstrongT=Astrong.';
LLSetup(lambda);

%% Make all unconnected points F points (if desired)
I=find(lambda==0);
II=[];
for i=I(:)'
  if isempty(find(AstrongT(:,i)))
    II(end+1)=i;
  end
end
I=II;
if length(I)>0
  if opt.make_unconnected_F
    disp(['   Making 2-way unconnected points F, number=',...
          num2str(length(I))])
    type(I)=F;
    LLRemove(I);
  else
    disp(['   Warning, found 2-way unconnected points, number=',...
          num2str(length(I))])
  end
end

while ~LLIsEmpty
  %% Make point with maximum lambda a coarse point
  i=LLGetMax;

  % Make this a C point
  type(i)=C;
  LLRemove(i);

  %% Make all remaining ST neighbours F points
  %STiU=find(Astrong(:,i) & type==U);
  STiU=CombiFind(Astrong(:,i), type, U);
  type(STiU)=F;
  LLRemove(STiU);

  %% Increase lambda by 1 for all remaining S neigbours of STiU
  for j=STiU(:)'
    %SjU=find(AstrongT(:,j) & type==U);
    SjU=CombiFind(AstrongT(:,j), type, U);
    %lambda(SjU)=lambda(SjU)+1;
    LLIncrease(SjU);
  end

  %% Decrease lambda by 1 for all remaining S neigbours of i
  %SiU=find(AstrongT(:,i) & type==U);
  SiU=CombiFind(AstrongT(:,i), type, U);
  %lambda(SiU)=lambda(SiU)-1;
  LLDecrease(SiU);
end
LLDestroy;

%%%% Cheat :

if isfield(opt,'perm') && length(opt.perm)==n
  disp('   Using given opt.perm and opt.nc')
  type(opt.perm(1:opt.nc))=C;
  type(opt.perm(opt.nc+1:n))=F;
end

%%%% A3 : Test for F-F connections without common C point

if opt.check_F_F>=1
  disp('   Checking for F-F connections without common C point')
  type=check_F_F(Astrong,AstrongT,type);
end
if opt.check_F_F>=2
  disp('   Checking for transpose F-F connections without common C point')
  type=check_F_F(AstrongT,Astrong,type);
end

%%%% Make perm and nc such that type(perm)=[C*ones(1,nc),F*ones(1,nc)]

nc=length(find(type==C));
[dummy,perm]=sort(type);

%%%% Error check

if length(find(type==U))>0
  disp('Error SelectCoarseGrid : U points left');
  type
end


%%%%%%%%%%%%%%%

function Keep=ThetaColumnCheck(As,opt)

n=size(As,2);

if opt.use_abs
  disp('   Using abs for Astrong')
  tol=opt.theta*max(abs(As));
  for i=1:n
    R=find(As(:,i));
    S=find(abs(As(R,i))>=tol(i));
    if length(S)>0
      As(R(S),i)=NaN;
    end
  end
else
  disp('   Using only negatives for Astrong')
  tol=opt.theta*full(max(-(As)));
  for i=1:n
    R=find(As(:,i));
    S=find(-(As(R,i))>=tol(i));
    if length(S)>0
      As(R(S),i)=NaN;
    end
  end
end

Keep=isnan(As);

return


%%%%%%%%%%%%

function I=CombiFind(A,B,c)

  %I=find(A & B==c);

  R=find(A);
  if length(R)>0
    I=R(B(R)==c);
  else
    I=[];
  end

return


%%%%%%%%%%%%

function type=check_F_F(Astrong,AstrongT,type)

  U=3;
  F=2;
  C=1;

  for i=find(type'==F)
    % still an F point ?
    if type(i)==F
      %Fi=find(AstrongT(:,i) & type==F);
      Fi=CombiFind(AstrongT(:,i), type, F);
      if length(Fi)>0
        %Ci=find(AstrongT(:,i) & type==C);
        Ci=CombiFind(AstrongT(:,i), type, C);
        problems=[];
        for j=Fi(:)'
          if length(find(Astrong(j,Ci)))==0
            % Problem
            problems(end+1)=j;
          end
        end
        if length(problems)>1
          % make i a C point
          type(i)=C;
        elseif length(problems)==1
          % make problem point a C point
          type(problems)=C;
        end
      end
    end
  end

return

%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%

% Linked lists for lambda's
%
% data structure
%
% LL.tbn(1:lambdamax,1:3) : [ top, bottom, number ]
% LL.lpn(1:3,1:n)         : [ lambda, previous, next ]
%
% We use increased lambda's to be able to have real negative lambda's

function LLSetup(lambda)
  global LL_tbn LL_lpn
  LL_tbn=[];
  LL_lpn=[];

  lambdamax=4*max(lambda);
  lambdamin=0;
  lambda=lambda-lambdamin+1;
  lambdamax=lambdamax-lambdamin+1;
  LL_tbn(1:lambdamax,1:3)=0;
  LL_lpn(1:3,1:length(lambda))=0;

  for i=1:length(lambda)
    LLAdd(i,lambda(i));
  end

return

function LLDestroy
  global LL_tbn LL_lpn

  clear global LL_tbn LL_lpn

return

function LLAdd(i,lambda)
  global LL_tbn LL_lpn

  if LL_tbn(lambda,3)==0
    % New
    LL_lpn(1:3,i)=[lambda;NaN;NaN];
    LL_tbn(lambda,1:2)=[i,i];
  else
    oldtop=LL_tbn(lambda,1);
    LL_lpn(:,i)=[lambda;NaN;oldtop];
    LL_lpn(2,oldtop)=i;
    LL_tbn(lambda,1)=i;
  end
  LL_tbn(lambda,3)=LL_tbn(lambda,3)+1;

return

function LLRemove(I)
  global LL_tbn LL_lpn

  for i=I(:)'
    lambda=LL_lpn(1,i);
    prev=LL_lpn(2,i);
    next=LL_lpn(3,i);

    if isnan(prev)
      LL_tbn(lambda,1)=next;
    else
      LL_lpn(3,prev)=next;
    end
    if isnan(next)
      LL_tbn(lambda,2)=prev;
    else
      LL_lpn(2,next)=prev;
    end
    LL_tbn(lambda,3)=LL_tbn(lambda,3)-1;
  end

return

function LLIncrease(I)
  global LL_tbn LL_lpn

  lambdamax=size(LL_tbn,1);

  for i=I(:)'
    lambda=LL_lpn(1,i)+1;
    if lambda>lambdamax
      fprintf('   Warning: lambda is running out of range, lambda=%d, i=%d\n',...
              lambda,i);
    else
      LLRemove(i);
      LLAdd(i,lambda);
    end
  end

return

function LLDecrease(I)
  global LL_tbn LL_lpn

  for i=I(:)'
    lambda=LL_lpn(1,i)-1;
    if lambda<1
      fprintf('   Warning: lambda is running out of range, lambda=%d, i=%d\n',...
              lambda,i);
    else
      LLRemove(i);
      LLAdd(i,lambda);
    end
  end

return

function i=LLGetMax
  global LL_tbn LL_lpn

  I=find(LL_tbn(:,3)>0);
  lmax=I(end);
  i=LL_tbn(lmax,1);

return

function val=LLIsEmpty
  global LL_tbn LL_lpn
%[length(find(LL_tbn(:,3)>0)),sum(LL_tbn(:,3))]
  val=length(find(LL_tbn(:,3)>0))==0;

return


