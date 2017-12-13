%Stencil2Matrix		Make matrix with given 2D stencil
%
%usage : A=Stencil2Matrix(S,m,n)
%
%S   = stencil matrix, for instance [  0 -1  0 ]
%                                   [ -1  4 -1 ]
%                                   [  0 -1  0 ]
%m,n = size of the 2D domain
%A   = The returned matrix, size m*n x m*n
%
%Note : At the boundary the first element outside the domain will be set
%       to zero, further elements will be a reflection of the real domain
%       In other words, for first order operators you get Dirichlet 
%       boundary condition, for second order also a Neuman boundary
%       condition.


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




function A=Stencil2Matrix(S,m,n)

si=(size(S,1)-1)/2;
sj=(size(S,2)-1)/2;

A=sparse(m*n,0);
for j=1:n
  for i=1:m
    R=sparse(m+2*si,n+2*sj);
    I=[i-si:i+si]+si;
    J=[j-sj:j+sj]+sj;
    R(I,J)=S;
    if si>1
      R(:,1+si:2*si-1)=R(:,1+si:2*si-1)+R(:,si-1:-1:1);
      R(:,m+2:m+si)=R(:,m+2:m+si)+R(:,m+2*si:-1:m+si+2);
    end
    if sj>1
      R(1+sj:2*sj-1,:)=R(1+sj:2*sj-1,:)+R(sj-1:-1:1,:);
      R(n+2:n+sj,:)=R(n+2:n+sj,:)+R(n+2*sj:-1:n+sj+2,:);
    end
    R=R([1:m]+si,[1:n]+sj);
    A(:,end+1)=R(:);
  end
end


A=A.';

  
