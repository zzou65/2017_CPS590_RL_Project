%AMGBlindS	Operator for the AMG coarse gird correction error operator
%
%The AMG data structure is passed in the gobal L_AMGBlindData, so this
%operator can be used in functions like jdqr.
%
%usage : x=AMGBlindS(b[,askdim])
%
%b      : vector to operate on
%x      : (1-P*inv(Ac)*R*A)b where the inv(Ac) is approximated with a
%         V-Cycle over the lower levels
%askdim : if this is set, x will return the number of unknowns (for jdqr)
%
%See also  AMGBlindErrorOperator AMGBlindS


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


function x=AMGBlindK(b,askdim)

global L_AMGBlindData
L=L_AMGBlindData;

if nargin>1
  x=size(L{1}.A,2);
  return
end


% Do x = b - P inv(Ac) R A b

r=L{1}.A*b;
rc=L{1}.R*r;
xc=AMGVcycle(L,rc,2);
x=L{1}.P*xc;
x=b-x;

