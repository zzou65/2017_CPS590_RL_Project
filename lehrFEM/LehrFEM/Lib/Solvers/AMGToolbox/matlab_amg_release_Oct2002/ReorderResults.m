%ReorderResults	Reoder AMG RunTests results
%
%usage : Rout=ReorderResults(R[[,Rorder],truncate])
%
%R      = results struct to reorder
%Rorder = results struct, Options array or index set to oder by,
%         default is taken from the LoopOptions function
%truncate : if 1, all options in R but not in Rorder are skipped,
%           if 0 they are placed at the end, default=0
%Rout   = The reordered and optionally truncated results struct
%
%Note : This will leave empty rows in the results where Rorder has
%       options not in R. This should not present problems
%
%See also  MergeResults ResultSortNames SaveResults


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


function R3=ReorderResults(R2,R1,truncate)

if length(fieldnames(R2))>7
  error('More fields than expected')
end

if nargin<3
  truncate=0;
end
if nargin<2 | length(R1)==0
  R1Opt=LoopOptions;
else
  if ~isstruct(R1)
    J=1:length(R2.Names);
    I=R1;
    R3.Conv=R2.Conv(I,J);
    R3.Prop=R2.Prop(I,J);
    R3.Times=R2.Times(I,J);
    R3.Errors=R2.Errors(I,J);
    R3.Names=R2.Names(J);
    R3.Options=R2.Options(I);
    R3.Opt=R2.Opt(I);
    return
  else 
    if isfield(R1,'Opt')
      R1Opt=R1.Opt;
    else
      R1Opt=R1;
    end
  end
end

  Default=AMGDefaultOptions;

  tmp(2)=R2.Opt(1);
  Empty=tmp(1);
  clear tmp

  for j=1:length(R1Opt)
    R1O(j)=mergestruct(Default,R1Opt(j));
  end

  I=[];
  I2=[];
  s=length(R1Opt);
  for i=1:length(R2.Opt)
    if comparestruct(R2.Opt(i), Empty)
    else
      R2O=mergestruct(Default,R2.Opt(i));
      found=0;
      for j=1:length(R1Opt)
        if comparestruct(R1O(j),R2O)
          fprintf('    Options R(%d) and Roder(%d) correspond \n',i,j)
          I(end+1)=j;
          I2(end+1)=i;
          found=1;
          break
        end
      end
      if ~found
        if truncate
          fprintf('New option  R(%d) is skipped\n',i);
        else
          s=s+1;
          I(end+1)=s;
          I2(end+1)=i;
          fprintf('New option  R(%d) becomes R(%d) \n',i,s);
        end
      end
    end
  end

  if 0 & isfield(R1,'Names')
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
	fprintf('New name R2(%d) = %s becomes R1(%d) \n',i,R2.Names{i},J(j));
	s=s+1;
	J(i)=s;
      end
    end
  else
    J=1:length(R2.Names);
  end

J2=1:length(R2.Names);

R3.Conv(1:max(I),1:max(J))=NaN;
R3.Conv(I,J)=R2.Conv(I2,J2);
R3.Prop(I,J)=R2.Prop(I2,J2);
R3.Times(I,J)=R2.Times(I2,J2);
R3.Errors(I,J)=R2.Errors(I2,J2);
R3.Names(J)=R2.Names(J2);
R3.Options(I)=R2.Options(I2);
R3.Opt(I)=R2.Opt(I2);

