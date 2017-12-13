%OptionsTable	Make a .tex and .ps table of an options struct array
%
%usage : [T,tit]=OptionsTable(Options[,AbreviateTitles[,FileName]])
%
%This will generate a table of the options in tex and ps format and
%show the ps file on screen using ghostview.
%
%Options : Options struct array, for instance Options or Opt in the
%          RunTests restuls or the results from LoopOptions
%AbreviateTitles : Whether to abreviate column titles (0 means full
%                  titles, 1 meand abreviated titles), default=0
%FileName : Filename base to use for the .tex and .ps file, default='table'
%T,tit : Table data that is passed to Table2Tex to make the tex and ps files
%
%See also  Table2Tex


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



function [T,tit]=OptionsTable(Options,AbreviateTitles,FileName)

if nargin<3
  FileName='table';
end
if nargin<2
  AbreviateTitles=0;
end

% modification nov 5 1999
header=input(' Enter Table Header in quotes %s:');


Options=Options';
n=length(Options);

tit={};
for nr=1:length(Options)
  O=Options(nr);
  name1=fieldnames(O);
  for i=1:length(name1)
    c=getfield(O,name1{i});
    if isstruct(c)
      name2=fieldnames(c);
      for ii=1:length(name2)
        cc=getfield(c,name2{ii});
        if isstruct(cc)
	  error('Third level of structs not allowed')
        else
          tit=AddToList(tit,[name1{i},'.',name2{ii}]);
        end
      end
    elseif length(c)>0 
      tit=AddToList(tit,[name1{i}]);
    end
  end
end

T={};
for i=1:length(Options)
  O=Options(i);
  for j=1:length(tit)
    T{i,j}=eval(['O.',tit{j}],'[]');
  end
end

%% Feed it to Table2Tex

% uses modified 
ifont=1;
ifont=input(' Enter 1 to use large font')
if ifont == 0
mTable2Tex(T,header,tit,AbreviateTitles,FileName);
elseif (ifont == 1)
bigfontmTable2Tex(T,header,tit,AbreviateTitles,FileName);
end;

return



%%%%%

function L=AddToList(L,name)

  if ~any(strcmp(L,name))
    L{end+1}=name;
  end
 
return
