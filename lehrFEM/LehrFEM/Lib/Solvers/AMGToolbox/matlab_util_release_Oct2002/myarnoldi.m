%myarnoldi	Very simple Arnoldi algorithem
%
%usage : [H,V]=myarnoldi(A,M,Mexplicit,V,nr);
%
%A : Matrix A or function for matvec of type y=A(x)
%V : initial vector
%nr: number of steps
%H : nr x nr projected matrix of S, H=V'*A*V
%V : Arnoldi basis, V'*V=I
%
%See also  arnoldi


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



function [H,V]=myarnoldi(A,V,nr);

k=0;
V=V/norm(V);
while k<nr
  % Arnoldi :
  % expand Krylov search (and schadow) space 
  % [V,H]=arnoldi(A,V,H);
  if ischar(A)
    v=eval([A,'(V(:,k+1))']);
  else
    v=A*V(:,k+1);
  end
  [V(:,k+2),H(1:k+2,k+1)]=mgs2(v,V);
  %[v,h]=mgs(V,v);
  %% h=V'*v;
  %% v=v-V*h;
  %hh=norm(v);
  %V=[V, v/hh];
  %H=[H, h ; zeros(1,k), hh ];
  k=k+1;
end

H=H(1:end-1,:);
V=V(:,1:end-1);

return

%%%%%%%%%%%

%MGS2		Modified Gram-Schidt orthonormalisation
%
%function [v,h]=mgs2(x,V,aant)
%
%input : V an n x j matrix with V'*V=I
%        x an n vector
%        aant : the number of times to repeat the othogonalisation
%output: v is x made otrhogonal to V, V'*v=0 and normalised
%        h is the j+1 vector containing the coefficients involved
%        x=[V,v]*h
%
%Especially usefull for the Arnoldi proces

function [v,h]=mgs2(x,V,aant)

% More stable version of GS :
% h=V'*v;
% v=v-V*s;

if nargin<3
  aant=1;
end

j=size(V,2);
h=zeros(j+1,1);
v=x;
for k=1:aant
  for i=1:j
    t=V(:,i)'*v;
    h(i)=h(i)+t;
    v=v-V(:,i)*t;
  end
end

h(j+1)=norm(v);
v=v/h(j+1);

%mgs2nul=norm(x-[V,v]*h)
%mgs2nul=norm(V'*v)
