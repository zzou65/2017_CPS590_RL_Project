%comparestruct	Test for equality of 2 structs
%
%usage val=comparestruct(a,b)
%
%a,b : The structs to compare
%val : value of a==b (1 for a==b, 0 for a~=b)


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



function val=comparestruct(a,b)

if ~isstruct(a) | ~isstruct(b)
  error('Input must be structs');
end

val=0;

aname=fieldnames(a);
bname=fieldnames(b);

if length(aname)~=length(bname)
  return
end

aname=sort(aname);
bname=sort(bname);

for i=1:length(aname)
  if ~strcmp(aname{i},bname{i})
    return
  end
end

for i=1:length(aname)
  av=getfield(a,aname{i});
  bv=getfield(b,bname{i});
  if xor( isstruct(av), isstruct(bv) )
    return
  end
  if isstruct(av)
    if ~comparestruct(av,bv)
      return
    end
  else
    if norm(size(av)~=size(bv),1)
      return
    end
    if length(av)>0 & norm(av~=bv,1)
      return
    end
  end
end

val=1;

return
