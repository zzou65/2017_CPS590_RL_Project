%arnoldi		Block Arnoldi recursion
%
%usage: [H,Q,b] = arnoldi(A0,vs,q,ijob);
%
%   Block Arnoldi recursion:  
%   qact = number of columns in Q on termination
%   On termination size of H is qact x qact - b
%   Q must contain at least q+b vectors before procedure terminates
%    where b may change during the computations but is the current
%    block size
%
%   Initial block size is size(vs,2).


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

  %%  Authors: Tong Zhang and Jane Cullum
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [H,Q,b] = arnoldi(A0,vs,q,ijob);

% Determine size of original system
% number of columns 
  nA=size(A0,2);

% track the block size in parameter b


% Check for deflation of starting block of vectors
%
 [u,t]=deflate(vs);

% Block size: Can change as we iterate
 b=size(u,2);

 H=zeros(q+b,q);
 Q=zeros(nA,q+b);


 s=b;
 pblk=1:b;
 Q(:,pblk)=u;
 bold = b;

  while s < q+b
% compute block or residue
% for H with q columns we need q+b Q-vectors
%
%  r=A*u if ijob=0; r=A-transpose if job=1
%  where A is original matrix. 
%
   if ijob == 0
     r = A0*u;
   else
     r=A0.'*u;
   end;
%

% reorthogonalize w.r.t previous q  
  for i=1:s
    H(i, pblk)=Q(:,i)'*r;
    r=r-Q(:,i)*H(i, pblk);
  end
%
%  generate next block

    [u,tt]=deflate(r);

    b=size(u,2);
    if b ~= bold
     bold = b;
    end;

    H(s+1:s+b,pblk)=tt;
    pblk=s+1:s+b;
    s=s+b;
    Q(:,pblk)=u;

% end of while 
 end
% need to return with qact which may be larger than q
%  depending upon the block sizes
qact=s;


% specify sizes of H and Q because I am looping on their sizes
H=H(1:s-b,1:s-b);
Q=Q(1:nA,1:s-b);
%











