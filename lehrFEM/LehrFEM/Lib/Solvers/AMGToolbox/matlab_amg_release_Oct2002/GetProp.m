%GetProp		Get property matrix from AMG RunTests Prop structure
%
%usage : p=GetProp(Prop,level,name)
%
%p(i,j)=Prop{i,j}(level).name
%
%if level=0 and name='levels' p(i,j)=number of levels for Prop(i,j)



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



function p=GetProp(P,l,name)

for i=1:size(P,1)
  for j=1:size(P,2)
    if strcmp( name, 'levels')
      if length(P{i,j})>0
        p(i,j)=length(P{i,j})+1;
      else
        p(i,j)=NaN;
      end
    else
      if length(P{i,j})>=l
        if ~iscell(P{i,j})
          t=eval(['P{i,j}(l).',name],'NaN');
        else
          t=eval(['P{i,j}{l}.',name],'NaN');
        end
        if length(t)>0
          p(i,j)=t(1);
        else
	  p(i,j)=NaN;
        end
      else
        p(i,j)=NaN;
      end
    end
  end
end
