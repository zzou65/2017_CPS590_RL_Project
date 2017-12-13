%changebasis	Change AMG data structure to P and R basis
%
%usage : L=changebasis(L)
%
%Changes everything to the PP=[P,[0;I]] and RR=[R;[0,I]] basis


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


function L=changebasis(L)

c=L{1}.c;
f=L{1}.f;
nc=length(c);
nf=length(f);
n=nc+nf;
cc=1:nc;
ff=nc+1:n;

PP(1:n,cc)=L{1}.P;
PP(f,ff)=speye(nf);
PP(c,ff)=0;

RR(cc,1:n)=L{1}.R;
RR(ff,f)=speye(nf);
RR(ff,c)=0;

L{1}.A=RR*L{1}.A*PP;
L{1}.P(cc,cc)=speye(nc);
L{1}.P(ff,cc)=0;
L{1}.R=L{1}.P';

L{1}.Mpre=inv(PP)*L{1}.Mpre*inv(RR);
L{1}.Mpost=inv(PP)*L{1}.Mpost*inv(RR);

L{1}.grid=L{1}.grid(L{1}.perm);
L{1}.c=cc;
L{1}.f=ff;
L{1}.perm=1:n;
L{1}.invperm=1:n;

L{1}.PP=PP;
L{1}.RR=RR;
