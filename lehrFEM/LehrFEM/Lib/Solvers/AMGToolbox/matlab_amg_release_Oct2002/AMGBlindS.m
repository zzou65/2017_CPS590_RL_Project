%AMGBlindS	Operator for the AMG smoother error operator Se
%
%The AMG data structure is passed in the gobal L_AMGBlindData, so this
%operator can be used in functions like jdqr.
%
%usage : x=AMGBlindS(b[,askdim])
%
%b      : vector to operate on
%x      : (1-Mpre*A)^its b
%askdim : if this is set, x will return the number of unknowns (for jdqr)
%
%See also  AMGBlindErrorOperator AMGBlindK


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



function x=AMGBlindS(b,askdim)

global L_AMGBlindData
L=L_AMGBlindData;

if nargin>1
  x=size(L{1}.A,2);
  return
end

x=iterate(L{1}.A,L{1}.Mpre,L{1}.Mpre_explicit,L{1}.A*b,L{1}.opt.pre.its);
x=b-x;

%%%%%%%%%%

function x=iterate(A,M,explicit,b,num)

if explicit
  x=M*b;
else
  x=M\b;
end
for i=2:num
  r=b-A*x;
  if explicit
    x=x+M*r;
  else
    x=x+M\r;
  end
end

return
