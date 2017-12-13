%MergeResults	Merge RunTests results
%
%usage : R=MergeResults(R1,R2)
%
%R1,R2 : Results from a RunTests run loaded in a struct with for instance
%        R1=load('Results.mat') or filenames of result files.
%        With only one argument this should either be a string with
%        wildcards indicating the filenames to use or a cell array of
%        filenames.
%R     : Struct of the same type with the combined results. Options and 
%        matrix names are matched, the order of R1 is retained for R.
%        Results in R2 overwrite results in R1.
%
%See also  ReorderResults ResultSortNames SaveResults


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


function R1=MergeResults(R1,R2)

if nargin<2
  if ischar(R1)
    d=dir(R1);
    R1=sort({d.name});
  end
  if ~iscell(R1)
    error('Argument error');
  end
  disp(['Starting with ',R1{1}]);
  R=load(R1{1});
  for i=2:length(R1)
    disp(['Merging ',R1{i}]);
    R=MergeResults(R,R1{i});
  end
  R1=R;
  return
else
  if ischar(R1)
    R1=load(R1);
  end
  if ischar(R2)
    R2=load(R2);
  end
end

  tmp(2)=R2.Opt(1);
  Empty=tmp(1);
  clear tmp

I1=[];
I2=[];
s=length(R1.Opt);
for i=1:length(R2.Opt)
  if comparestruct(R2.Opt(i), Empty)
    fprintf('Option  R2(%d) is empty \n',i);
  else
    found=0;
    for j=1:length(R1.Opt)
      if comparestruct(R1.Opt(j),R2.Opt(i))
        fprintf('Options R1(%d) and R2(%d) correspond \n',j,i)
        I1(end+1)=j;
        I2(end+1)=i;
        found=1;
        break
      end
    end
    if ~found
      s=s+1;
      I1(end+1)=s;
      I2(end+1)=i;
      fprintf('New option R2(%d) becomes R1(%d) \n',i,I(i));
    end
  end
end

s=length(R1.Names);
for i=1:length(R2.Names)
  found=0;
  for j=1:length(R1.Names)
    if strcmp(R1.Names{j},R2.Names{i})
      fprintf('Names R1(%d) and R2(%d) correspond = %s \n',j,i,R2.Names{i});
      J(i)=j;
      found=1;
      break
    end
  end
  if ~found
    s=s+1;
    J(i)=s;
    fprintf('New name R2(%d) = %s becomes R1(%d) \n',i,R2.Names{i},J(i));
  end
end

oldsize1=size(R1.Conv,1);
oldsize2=size(R1.Conv,2);

for i=1:length(I1)
  for j=1:size(R2.Conv,2)
    if length(R2.Prop{I2(i),j})>0 | ischar(R2.Errors{I2(i),j})
      if I1(i)<=oldsize1 & J(j)<=oldsize2 & length(R1.Prop{I1(i),J(j)})>0
        fprintf('Warning : Overwriting R1(%d,%d) with R2(%d,%d)\n',...
                 I1(i),J(j),I2(i),j);
      else
        fprintf('              Writing R1(%d,%d) with R2(%d,%d)\n',...
                 I1(i),J(j),I2(i),j);
      end
      R1.Conv(I1(i),J(j))=R2.Conv(I2(i),j);
      R1.Prop{I1(i),J(j)}=R2.Prop{I2(i),j};
      R1.Times(I1(i),J(j))=R2.Times(I2(i),j);
      R1.Errors{I1(i),J(j)}=R2.Errors{I2(i),j};
    end
  end
end
for j=1:size(R2.Conv,2)
  R1.Names{J(j)}=R2.Names{j};
end
for i=length(I1)
  R1.Options(I1(i))=SortStruct(R2.Options(I2(i)),fieldnames(R1.Options(1)));
  R1.Opt(I1(i))=R2.Opt(I2(i));
end


%%%

function s=SortStruct(sin,names)

s=[];
for i=1:length(names)
  s=setfield(s,names{i},getfield(sin,names{i}));
end
