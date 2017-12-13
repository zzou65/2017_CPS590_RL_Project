%sortbases	Sort a base and possibly make vectors as real as possible
%
%function [e,v]=sortbase(e,v,operation,rotate,normalise)
%
% operation(e) is the vector of values to sort by
% v is the matrix (set of columns) to be sorted
% optional : rotate>=1 results in "realifying" the vectors
% optional : normalise>=1 results in normalisation of the vectors


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



function [e,v]=sorteig(e,v,operation,rotate,normalise)

if (nargin<5)
  normalise=0;
end
if (nargin<4)
  rotate=0;
end
if nargin<3
  operation='abs';
end

[dummy,I]=eval([ 'sort(',operation,'(e))' ]);
e=e(I);
v=v(:,I);

if (rotate>0)
  I=sqrt(-1);
  for i=1:size(v,2)
    t=v(find(abs(v(:,i))>0),i);
    v(:,i)=v(:,i).*exp(-I*sum(atan(imag(t)./real(t)).*abs(t))/sum(abs(t)));
    if (real(v(1,i))<0)
      v(:,i)=-v(:,i);
    end
  end
end

if (normalise>0)
  for i=1:size(v,2)
    v(:,i)=v(:,i)./norm(v(:,i));
  end
end


