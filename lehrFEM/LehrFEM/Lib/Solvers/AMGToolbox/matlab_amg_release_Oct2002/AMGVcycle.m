%AMGVcycle	Preform an AMG V-cycle
%
%   X = AMGVCYCLE(AMGDATA,B) performs algebraic mutlgird V-cycle to
%   approxmate solution X of
%       AMGDATA{1}.A*X = B
%   starting from a zero initial guess.  To use an arbitrary initial guess
%   X0, use (for A = AMGDATA{1}.A)
%       C = AMGVCYCLE(AMGDATA, B - A*X0);
%       X = X + C;
%   AMGDATA is the AMG data structure generated with AMGSetup.
%
%   Original description: -----------------------------------------------|
%   |                                                                    |
%   | usage : x=AMGVcycle(L,b[,level])                                   |
%   |                                                                    |
%   | L     = AMG data structure from AMGSetup                           |
%   | b     = residual to solve for                                      |
%   | level = level we are at, default=1, for recursive use only         |
%   | x     = approximate solution to L{l}.A x = b                       |
%   |____________________________________________________________________|
%
% 	See also  AMGSetup readall amg_solve namg_solve


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

function x=AMGVcycle(L,b,l)

% Initialize

global ndffiters cmaxiters;
%

idebug=0;
if idebug == 1;  disp(' Entering AMGVcycle: level ='); disp(l);end;

opt=L{1}.opt;

if nargin<3
  l=1;
end

if l<length(L)
  %% Regular level

  if idebug==1; disp('Regular Level: before Presmooth'); end;
  if idebug==1; disp('opt.pre.its');disp(opt.pre.its);end;

  % Pre-Smooth
  if opt.pre.its == 0
    % No Pre-Smoothing
    nx=length(b);
    x1=zeros(nx,1);
  else
    % Apply specified smoother (with zero initial guess)
    x1=iterate(L{l}.A,L{l}.Mpre,L{l}.Mpre_explicit,b,opt.pre.its);
  end;

  %if idebug==1; disp('before r1 computation');disp('b');disp(b(1:10));end;
  if idebug==1; disp('before r1 computation');disp('b');disp(b);end;

  r1=b-L{l}.A*x1; % residual on fine mesh

  %if idebug== 1;disp('r1');disp(r1(1:10));disp(' before restriction of r');end;
  if idebug== 1;disp('r1');disp(r1);disp(' before restriction of r');end;

  % Coarse grid correction
  rc=L{l}.R*r1;   % residual restricted to coarse mesh

  %if idebug==1; disp('rc');disp(rc(1:10));end;
  if idebug==1; disp('rc');disp(rc);end;
  if idebug==1; disp('before recursive call to AMGVcycle: level l+1');disp(l+1);end;

  xc=AMGVcycle(L,rc,l+1); % coarse grid correction

  if idebug==1; disp('before x2 computation');disp('xc');disp(xc);end;
  %if idebug==1; disp('before x2 computation');disp('xc');disp(xc(1:10));end;

  x2=x1+L{l}.P*xc;  % prolong correction to fine grid and update solution

  % jane modification nov 17 1999
  % below replaces the brute force 1st attempt in nov 11 version
  % the alternative would be to set mDff to zero and not do test below.
  if L{l}.imDff ~= 0
    %    x2(L{l}.f) = x2(L{l}.f) + L{l}.mDff \ r(L{l}.f);
    % for now we leave as subroutine so we can also do exact correction.
    %  i.e. if opt.method == 'exact' and L{l}.imDff== 1 then use exact correction
    x2=ModifyDelta(L,r1,l,x2);
  end;

  %%

  if idebug==1; disp('before r2 computation');disp('x2');disp(x2);end;
  %if idebug==1; disp('before r2 computation');disp('x2');disp(x2(1:10));end;

  r2=b-L{l}.A*x2; % residual : right-hand side in postsmoothing
  % Post-Smooth
  %if idebug==1;disp('r2');disp(r2(1:10));end;
  if idebug==1;disp('r2');disp(r2);end;
  if idebug==1; disp(' before post smoothing: Mpost: Mpost_explicit');end;
  if idebug==1; disp('opt.post.its');disp(opt.post.its);end;

  if opt.post.its == 0
    % no post smoothing
    x=x2;
  else % do postsmoothing with residual as right hand side and zero initial guess
    x=x2+iterate(L{l}.A,L{l}.Mpost,L{l}.Mpost_explicit,r2,opt.post.its);
  end;
  %if idebug==1;disp('x');disp(x(1:10));end;
  if idebug==1;disp('x');disp(x);end;
  if idebug==1; disp('size of x');disp(size(x));end;

  %x=[rc,xc];

else % coarsest level

  if idebug==1; disp(' at bottom level: explicit solve');end;
  if idebug==1; disp(' level'); disp(l);end;

  % Bottom level
  % solve equation exactly using (precomputed) LU-factorization
  x=L{l}.U\(L{l}.L\b);
  %x=L{l}.A\b;
  %if idebug==1; disp('x');disp(x(1:10));end;
  if idebug==1; disp('x');disp(x);end;
end

if idebug==1; disp(' at end of AMGVcylce subroutine');end;

return

%%%%%%%%%% smoothing iteration

function x=iterate(A,M,explicit,b,num)
%
% A is the stiffness matrix, M is an approximate inverse of A if explicit
% is true and an approximation of A if explicit is false, b is the right
% hand side and num is the number of iterations.

idebug=0;
if idebug==1;disp(' in subroutine iterate');end;

if explicit

  if idebug==1; disp(' compute x=M*b in iterate: M*A*e');end;

  x=M*b;
else

  if idebug==1; disp(' compute M\b in smoother: invQ*A*e');end;

  x=M\b;
end

for i=2:num

  if idebug==1; disp(' compute new residual in iterate:');end;
  if idebug==1; disp(' r=(I - A*inv(Q))r: e=(I-inv(Q)*A)');end;

  r=b-A*x;

  if explicit

    if idebug==1; disp(' in explicit in i loop in iterate: x=x+M*A*e') ;end;

    x=x+M*r;
  else

    if idebug==1; disp(' in not explicit in i loop in iterate: x=x+inv(Q)*r');end;

    x=x+M\r;
  end
end

if idebug==1; disp(' iterate call finished');end;
return

%%%%%%%%%%% extension to AMG : extra correction of fine points after coarse
%%%%%%%%%%% grid correction.

function [modx]=ModifyDelta(L,r,l,x)

% Adds correction to basic AMG Correction x: Exact correction would
%  be A_{FF}^{-1}*r_F term. we are not trying to
%  approximate it but rather to help the smoother by
%  forcing corrections to r_F in addition to those made by the smoother.

%  modx is x + (correction on fine points only).

%
% pick off the fine points
f=L{l}.f;
fr = r(f);

% only enter this subroutine if opt.P.mDff = 1 so
%  we are supposed to do somesort of correction

if L{l}.exact ~= -1
  % Use L{l}.mDff
  Dff = L{l}.mDff;
  disp([' Using Dff to Correct: level=',num2str(l)]);
else     % L{l}.exact = -1
  % Use exact Aff
  Dff = L{l}.A(f,f);
  disp([' Exact: Using Aff to Correct: level=',num2str(l)]);
end;
%
fdelta= Dff\fr;
nfdelta=norm(fdelta);
disp(' norm of correction');disp(nfdelta);
%
modx = x;
modx(f) = modx(f) + fdelta;
%
return
%%%
