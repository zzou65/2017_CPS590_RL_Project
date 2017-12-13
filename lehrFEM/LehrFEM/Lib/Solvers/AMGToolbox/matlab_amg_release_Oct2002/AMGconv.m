%AMGconv		Get maximum eigenvalue of overall AMG erorr operator
%
%usage : c = AMGconv(L[,l])
%
%L     : AMG data structure from AMGSetup
%level : Level to use, default=1
%c     : maximum eigenvalue of overall AMG erorr operator, cumputed
%        using JDQR, if jdqr does not converge [] is returned
%
%See also  AMGBlindErrorOperator jdqr


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



function r = AMGconv( L, l )

if nargin<2
  l=1;
end

global L_AMGBlindData
L_AMGBlindData={L{l:end}};
L_AMGBlindData{1}.opt=L{1}.opt;

%[ev,e,h]=jdqr('AMGBlindErrorOperator',1);
%r=max(diag(e));
%plot(h(:,3),log10(h(:,1)),'-o')

r=jdqr('AMGBlindErrorOperator',1);

clear global L_AMGBlindData


