%gsai		Compute generalised sparse approximate inverse
%
%function [M,res]=gsai(A,B,S)
%
%Make a sparse approx. to inv(A)*B, min || AM-B ||F
%
%A,B		input matrices, A must be square ??
%S		matrix which sparsity patern is the sp.pat. for M
%
%M		approximate inverse of A
%res		Vector of column residuals



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



function [M,res]=gsai(A,B,S)

if size(A,1)~=size(B,1)
  error('size mismatch')
end

l=size(A,1);
m=size(A,2);
n=size(B,2);
M=spconvert([m,n,0]);	% Empty nxn matrix

if ( nnz(A)/(n*n) < 0.30 )
  Afull=0;
  Sa=spones(A);
  %Sb=spones(B);
else
  Afull=1
end

Tindex=0;
Tsolve=0;


for k=1:n

  % Get sparsity patern for this column
  J=find(S(:,k))';

  tic;
  I=[];
  if Afull
    I=[1:n];
    %Ib=I;
  else
    SaSum=sum(Sa(:,J),2);
    I=find(SaSum);
    %Ib=find(SaSum+Sb(:,k));
  end
  Tindex=Tindex+toc;

  %disp([ k, length(I), length(J) ])

  tic;
  %r=B(:,k);
  if length(J)>0
%    M(J,k)=full(A(I,J))\B(I,k);
    M(J,k)=gelss(full(A(I,J)),B(I,k),1e-10);
    %r(Ib)=r(Ib)+A(Ib,J)*M(J,k);
  end
  %res(k)=norm(r);
  Tsolve=Tsolve+toc;
end
       
Tindex
Tsolve



%GELSS                 Solve LS problem for possibly rank deficient A
%
%function [x,rank,res]=gelss(A,b,rcond)

function [x,rank,res]=gelss(A,b,rcond)

if length(A)>0

if nargin<3
  rcond=1e-12;
end
residu=(nargout>2);

m=size(A,1);
n=size(A,2);

[U,S,V]=svd(A,0);
if min(m,n)>1
  S=diag(S);
end
h=find( S >= S(1)*rcond );
rank=length(h);

d=U'*b;
y=d(h)./S(h);
x=V(:,h)*y;

if residu
  res=norm(A*x-b);
end

else

  x=zeros(size(A,2),size(b,2));
  rank=0;
  res=norm(b);

end


