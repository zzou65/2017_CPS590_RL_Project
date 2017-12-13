%ComputeProperties	Compute the properties of an AMG setup
% Modified Oct/Nov 1999
%usage : P=ComputeProperties(L[,DoEig])
%
%L = AMG data struct form AMGSetup
%DoEig = if 1, does eigenvalue convergence analysis, if 0 not. Default
%        is to do it when the matrix is smaller than 5000
%P = Properties struct


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


function P=ComputeProperties(L,DoEig)

global L_AMGBlindData

if nargin<2 | length(DoEig)==0
  DoEig=(length(L{1}.A)<5000);
end
if DoEig
  disp('   Will do convergence eigenvalue analysis');
else
  disp('   Will not do convergence eigenvalue analysis');
end

P=[];
for l=1:length(L)-1
  A=L{l}.A;
  As=L{l}.Astrong;
  D=spdiag(diag(A));
  Ds=spdiag(diag(As));
  f=L{l}.f;
  c=L{l}.c;
  n=size(A,2);
  L_AMGBlindData={L{l:end}};
  L_AMGBlindData{1}.opt=L{1}.opt;


  if DoEig
    P(l).Convergence        = jdqr('AMGBlindErrorOperator',1);
    if length(P(l).Convergence)==0
      P(l).Convergence=NaN;
    end
  end
%   disp([' Before SolveConv:Level',num2str(l)]); 
%  [conv,itercpu,iterflops]= SolveConvergence(L,l);
%   disp([' After SolveConv:Level',num2str(l)]); 
% original is above
  [resid,conv,itercpu,iterflops]= SolveConvergence(L,l);
  P(l).SCresid = resid;


idebug=0;
if idebug == 1
 disp('in computeproperties after solveconvergence ')
end;
  P(l).WayneConvergence   = conv;
  P(l).IterationTime      = [itercpu,iterflops];
  P(l).TimePerWayneOrderPerN  = [itercpu,iterflops]*...
                                log(0.1)/log(P(l).WayneConvergence(1))/n;
  if DoEig
    P(l).TimePerOrderPerN   = [itercpu,iterflops]*...
                              log(0.1)/log(abs(P(l).Convergence))/n;
  end
  P(l).n               = length(A);
  P(l).nc              = length(c);
  P(l).nf              = length(f);
  P(l).A.NnzPerRow     = nnz(A)/size(A,1);
  P(l).As.NnzPerRow     = nnz(As)/size(A,1);
  P(l).A.Nnz     = nnz(A);
  P(l).As.Nnz     = nnz(As);
  P(l).A.AntiSymmetry  = AntiSymmetry(A);
  P(l).As.AntiSymmetry  = AntiSymmetry(As);
%  P(l).A.RangeRealEigs = RangeRealEigs(A);
  P(l).A.RangeRowsums  = Range(sum(A.'));
  P(l).As.RangeRowsums  = Range(sum(As.'));
  P(l).A.RangeDiag     = Range(diag(A));
  P(l).As.RangeDiag     = Range(diag(As));
  P(l).A.RangeOffDiag  = Range(A-D);
  P(l).As.RangeOffDiag  = Range(As-Ds);
%  P(l).A.IsMMatrix     = real(P(l).A.RangeRealEigs(1))>0 & ...
%                         P(l).A.RangeDiag(1)>0 & ...
%                         P(l).A.RangeOffDiag(2)<0 ;
if idebug ==1
  disp(' in computeproperties after P.as.rangeofdiag'); end;
 
% below menno original
%  P(l).A.RangeDelta    = Range( 1 - sum(abs(A-D)')./diag(D)' );
%  P(l).As.RangeDelta    = Range( 1 - sum(abs(As-Ds)')./diag(D)' );
% below jane modified
% cannot use Ds because with lumping some Ds can be 0
  P(l).A.RangeDelta    = Range( 1 - sum(abs(A-D)')./diag(abs(D))' );
  P(l).As.RangeDelta    = Range( 1 - sum(abs(As-Ds)')./diag(abs(D))' );
  P(l).Aff.RangeDiag   = Range(diag(A(f,f)));
  P(l).Aff.RangeOffDiag= Range(A(f,f)-D(f,f));
  P(l).Aff.RangeDelta    = Range( 1 - sum(abs(A(f,f) -D(f,f))')./diag(abs(D(f,f)))' );
  P(l).Aff.NnzPerRow   = nnz(A(f,f))/length(f);
  P(l).Asff.RangeDiag   = Range(diag(As(f,f)));
  P(l).Asff.RangeOffDiag= Range(As(f,f)-Ds(f,f));
  P(l).Asff.RangeDelta    = Range( 1 - sum(abs(As(f,f) -Ds(f,f))')./diag(abs(D(f,f)))' );
  P(l).Asff.NnzPerRow   = nnz(As(f,f))/length(f);
  P(l).Afc.RangeOffDiag= Range(A(f,c));
  P(l).Afc.NnzPerRow   = nnz(A(f,c))/length(f);
  P(l).Asfc.RangeOffDiag= Range(As(f,c));
  P(l).Asfc.NnzPerRow   = nnz(As(f,c))/length(f);
%  P(l).S.RangeRealEigs = RangeRealEigs('AMGBlindS');
% modification nov 2 1999
%iskip=0;
%if iskip ~= 0
%  [P(l).As.RangeFractionPosElements]=ComputeFraction(As,Ds);
%  [P(l).A.RangeFractionPosElements]=ComputeFraction(A,D);
%end;

end

%%%%%%%

function val=AntiSymmetry(A)
  val=max(max(abs(A-A.')))./max(max(abs(A+A')));
return

%%%%%%

function val=RangeRealEigs(A)
  %opt.MaxIt=25;
  %t=jdqr(A,1,'SR',opt);
  t=jdqr(A,1,'SR');
  if length(t)==0
    t=NaN;
  end
  val(1)=t;
  %t=jdqr(A,1,'LR',opt);
  t=jdqr(A,1,'LR');
  if length(t)==0
    t=NaN;
  end
  val(2)=t;
return

%%%%%%

function val=RangeAbsEigs(A)
  %opt.MaxIt=25;
  %t=jdqr(A,1,'SM',opt);
  t=jdqr(A,1,'SM');
  if length(t)==0
    t=NaN;
  end
  val(1)=t;
  %t=jdqr(A,1,'LM',opt);
  t=jdqr(A,1,'LM');
  if length(t)==0
    t=NaN;
  end
  val(2)=t;
return

%%%%%%

function val=Range(x)
  val(1)=min(x(:));
  val(2)=max(x(:));
return

%%%%%%
%jane: nov 19, 1999
function [resnorm,conv,itercpu,iterflops]=SolveConvergence(L,l)
% below is original
%function [conv,itercpu,iterflops]=SolveConvergence(L,l)


%  maxit=25;
  maxit=50;
  reltol=1e-7;
%  reltol=1e-10;
% Choose right-hand side
%  istart=0;
  istart=0;
  if istart == 0
   disp(' Using Random Right Hand Side');
% original menno
%   rand('state',130372);
   rand('state',7463193);
   b=rand(size(L{l}.A,2),1);
  elseif(istart == -1 & l==1 )
% read in b from file
   disp(' Reading in b from file');
   eval(['load rhs.dat']); 
   b=rhs;
  elseif(istart == 1 | (istart == -1 & l > 1))
   disp(' Using Random Xtrue: Compute b=A*Xtrue');
% modify to track errors and residuals
%   rand('state',130372);
   rand('state',7463193);
% generate xtrue  randomly
   xtrue=rand(size(L{l}.A,2),1);
   normxtrue=norm(xtrue);
%  xtrue=ones(size(L{l}.A,2),1);
% compute b
   b= L{l}.A*xtrue; 
  end;
idebug=0; if idebug==2; disp(xtrue);disp(size(L{l}.A,2));disp(b);end;
  normb=norm(b);
% track resnorm
  iresidual=0;
  resnorm=[];
%
idebug=0;
if idebug==1;disp('entering SolveConvergence: level');disp(l);end;
if idebug == 1;  disp(normb);end;
  x=zeros(size(b));
  r=b;
  k=0;
  normr=normb;
% added line
  resnorm=[resnorm normr/normb];
%
  startcpu=cputime;
  startflops=0; %flops;
% want to track the error on each level
%
  while normr>reltol*normb & k<maxit & normr<1e100
%
%% added and then deleted
%   ndffiters=k;
if idebug==1;xold=x;end;
    x=x+AMGVcycle(L,r,l);
if idebug==1; disp('x minus xold');test=max(abs(x-xold));disp(test);end;
    r=b-L{l}.A*x;
if idebug==1;disp(' r minus b');test2=max(abs(r-b));disp(test2);end;
    k=k+1;
    normrp=normr;
    normr=norm(r);
% modification nov 2 1999
    resnorm=[resnorm normr/normb];
%
    convk=normr/normrp;
if idebug==1;disp(k);disp(normr);disp(convk);end;
    if k==1
      conv1=convk;
    end
  end
  itercpu=(cputime-startcpu)/k;
  iterflops= 0; %(flops-startflops)/k;
  avrconv=(normr/normb)^(1/k);
  conv=[avrconv,conv1,convk];
if idebug ==1; disp(' at end of solveconvergence before plot');end;



if iresidual == 1;
figure;
plot(log10(resnorm),'b.');
title(' Residual Norms: Solve Convergence: ');
xlabel([' Iteration Number: Level:',num2str(l)]);
ylabel(' Log10(abs(||r(k)||/||r(0)||) ');
if idebug==1;disp(' end of solveconvergence after plot');end;
end;



return
