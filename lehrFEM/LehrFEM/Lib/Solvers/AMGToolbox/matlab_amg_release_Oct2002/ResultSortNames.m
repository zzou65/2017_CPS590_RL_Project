%ResultSortNames	Reoder AMG RunTests results
%
%usage : Rout=ResultSortNames(R[,J])
%
%R      = results struct to reorder
%J      = Order vector, default reverse alphabetic with increasing numbers
%Rout   = The reordered and optionally truncated results struct
%
%See also  MergeResults ReorderResults SaveResults


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


function Rout=ReorderResults(R,J)

if length(fieldnames(R))>7
  error('More fields than expected')
end

if nargin<2
  for i=1:length(R.Names)
    N{i}=FillNumbers(R.Names{i});
  end
  [dummy,J]=sort(N);
  J=J(end:-1:1);
end

I=1:length(R.Options);

Rout.Conv=R.Conv(I,J);
Rout.Prop=R.Prop(I,J);
Rout.Times=R.Times(I,J);
Rout.Errors=R.Errors(I,J);
Rout.Names=R.Names(J);
Rout.Options=R.Options(I);
Rout.Opt=R.Opt(I);

return



function t=FillNumbers(s)

  max=10;
  s=[s,' '];
  numbers=(s>='0')&(s<='9');
  start=find(numbers(1:end-1)==0 & numbers(2:end)==1)+1;
  stop =find(numbers(1:end-1)==1 & numbers(2:end)==0);
  last=0;
  t=[];
  for i=1:length(start)
    ln=stop(i)-start(i)+1;
    %ns=[char('0'*ones(1,max-ln)),s(start(i):stop(i))];
    ns=num2str(10^max-str2num(s(start(i):stop(i)))-1,max);
    t=[t,s(last+1:start(i)-1),ns];
    last=stop(i);
  end
  t=[t,s(last:end-1)];

return      
