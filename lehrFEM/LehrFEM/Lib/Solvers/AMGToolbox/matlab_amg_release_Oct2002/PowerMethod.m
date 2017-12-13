

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



function v=MPowerMethod(A,v,nr)

  for i=1:nr
% ischar(S) returns 1 if S is character array: 0 Otherwise
    if ischar(A)
% 
      v=eval([A,'(v)']);
    else
      v=A*v;
    end
    v=v/norm(v);
  end

return

%%%%%%%%%%%%


function v=MArnoldiMethod(A,v,nr)

  % Get Arnoldi AV~VH
  % Note: I get the square H
  %[H,V]=arnoldi(A,v,nr,0);
  [H,V]=myarnoldi(A,v,nr);
% Computes eigenvalues and eigenvectors of H
% and sorts them by abs(eval)
  [e,ev]=sorteig(H);
  v=V*ev(:,end);
% this assumes we have real eigenvectors
  v=real(v);

return
