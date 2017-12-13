%sorteig		Calculate eigenvalues and -vectors and sort them
%
%usage : [e,v]=sorteig(a[[,operation],rotate])
%
%This uses the matlab builtin function eig to get the eigenvalues and -vectors
%
%A         = Matrix to get eigenvalues and -vectors of
%operation = string with the function to use for the sort, default='abs'
%rotate    = if 1, the vectors will be multiplied with e^i\phi in an attempt
%            to maximize the real part of the vector, default=1
%e         = the vector of eigenvalues sorted by operation(e)
%v         = the matrix of corresponding eigenvectors, v=[v1,v2,..,vn]


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



function [e,v]=sorteig(a,operation,rotate)

if (nargin<2)
  operation='abs';
end
if (nargin<3)
  rotate=1;
end



a=full(a);
[v,e]=eig(a);
e=diag(e);
[e,v]=sortbase(e,v,operation,rotate);


