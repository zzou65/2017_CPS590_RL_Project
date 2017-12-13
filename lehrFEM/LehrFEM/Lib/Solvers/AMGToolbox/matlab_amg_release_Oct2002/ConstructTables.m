%ConstructTables   Make a .tex and .ps table of a matrix
%
%usage : [T,tit,M]=ConstructTables(Options,Prop,Names[,AbreviateTitles[,FileName]])
%
%This will generate a table in tex and ps format and
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


function [T,tit,M]=ConstructTables(Options,Prop,Names,AbreviateTitles,FileName)

if nargin<5
  FileName='table';
end
if nargin<4
  AbreviateTitles=0;
end

Options=Options';
no=length(Options);

clear M;
TMatNo=length(Names);

% property must be a string
% example 'A.RangeDiag'


inp=0;
while inp~=1&inp~=2&inp~=3&inp~=4&inp~=5&inp~=6&inp~=7&inp~=8&inp~=9&inp~=10&inp~=11&inp~=12&inp~=13&inp~=14&inp~=15&inp~=16&inp~=17&inp~=18&inp~=19&inp~=20&inp~=21&inp~=22&inp~=23&inp~=24&inp~=25&inp~=26&inp~=27&inp~=28&inp~=29&inp~=30&inp~=31&inp~=32&inp~=50&inp~=35&inp~=36;
  disp('Please select from the following properties');
  disp('1: A.NnzPerRow');
  disp('2: A.AntiSymmetry');
  disp('3: As.NnzPerRow');
  disp('4: As.AntiSymmetry');
  disp('5: Afc.NnzPerRow');
  disp('6: Asfc.NnzPerRow');
  disp('7: Aff.NnzPerRow');
  disp('8: Asff.NnzPerRow');
  disp('9: nc');
  disp('10: nf');
  disp('11: n');
  disp('12: A.RangeDiag');
  disp('13: A.RangeRowSums');
  disp('14: A.RangeOffDiag');
  disp('15: A.RangeFractionPosElements');
  disp('16: A.RangeDelta');
  disp('17: As.RangeDiag');
  disp('18: As.RangeRowSums');
  disp('19: As.RangeOffDiag');
  disp('20: As.RangeFractionPosElements');
  disp('21: As.RangeDelta');
  disp('22: Aff.RangeDiag');
  disp('23: Asff.RangeDiag');
  disp('24: Aff.RangeOffDiag');
  disp('25: Asff.RangeOffDiag');
  disp('26: Afc.RangeOffDiag');
  disp('27: Asfc.RangeOffDiag');
  disp('28: IterationTime');
  disp('29: TimePerWayneOrderPerN');
  disp('30: WayneConvergence: All 3: av:1st:last');
  disp('35: A.Nnz');
  disp('36: As.Nnz');
  disp('31: Aff.RangeDelta');
  disp('32: Asff.RangeDelta');
  disp('50: WayneConvergenceMatrix');
inp = input('Enter the number of the desired property matrix:');
end;

clear M;


for m=1:TMatNo;
  MatNo=m;

for k=1:no
% number of levels for option k
nl=size(Prop{k,MatNo},2);
 for j=1:nl;
  if inp == 1
  M(k,j,m)=Prop{k,MatNo}(j).A.NnzPerRow;
elseif inp == 2
  M(k,j,m)=full(Prop{k,MatNo}(j).A.AntiSymmetry);
elseif inp == 3
  M(k,j,m)=Prop{k,MatNo}(j).As.NnzPerRow;
elseif inp == 4
  M(k,j,m)=full(Prop{k,MatNo}(j).As.AntiSymmetry);
elseif inp == 5
  M(k,j,m)=Prop{k,MatNo}(j).Afc.NnzPerRow;
elseif inp == 6
  M(k,j,m)=Prop{k,MatNo}(j).Asfc.NnzPerRow;
elseif inp == 7
  M(k,j,m)=Prop{k,MatNo}(j).Aff.NnzPerRow;
elseif inp == 8
  M(k,j,m)=Prop{k,MatNo}(j).Asff.NnzPerRow;
elseif inp == 9
  M(k,j,m)=Prop{k,MatNo}(j).nc;
elseif inp == 10
  M(k,j,m)=Prop{k,MatNo}(j).nf;
elseif inp == 11
  M(k,j,m)=Prop{k,MatNo}(j).n;
elseif inp == 50
  M(k,j,m)=Prop{k,MatNo}(j).WayneConvergence(1);
elseif inp == 35
  M(k,j,m)=Prop{k,MatNo}(j).A.Nnz;
elseif inp == 36
  M(k,j,m)=Prop{k,MatNo}(j).As.Nnz;
 end;
end;
end;


for k=1:no
% number of levels for option k
 nl=size(Prop{k,MatNo},2);
 i1=0; i2=0;
 for j=1:nl;
  i1=i2+1;i2=i1+1;
  if inp == 12
   M(k,i1:i2,m)=full(Prop{k,MatNo}(j).A.RangeDiag);
  elseif inp == 13
   M(k,i1:i2,m)=full(Prop{k,MatNo}(j).A.RangeRowsums);
  elseif inp == 14
  M(k,i1:i2,m)=full(Prop{k,MatNo}(j).A.RangeOffDiag);
elseif inp == 15
  M(k,i1:i2,m)=full(Prop{k,MatNo}(j).A.RangeFractionPosElements);
elseif inp == 16
  M(k,i1:i2,m)=full(Prop{k,MatNo}(j).A.RangeDelta);
elseif inp == 17
  M(k,i1:i2,m)=full(Prop{k,MatNo}(j).As.RangeDiag);
elseif inp == 18
  M(k,i1:i2,m)=full(Prop{k,MatNo}(j).As.RangeRowsums);
elseif inp == 19
  M(k,i1:i2,m)=full(Prop{k,MatNo}(j).As.RangeOffDiag);
elseif inp == 20
  M(k,i1:i2,m)=full(Prop{k,MatNo}(j).As.RangeFractionPosElements);
elseif inp == 21
  M(k,i1:i2,m)=full(Prop{k,MatNo}(j).As.RangeDelta);
elseif inp == 22
  M(k,i1:i2,m)=full(Prop{k,MatNo}(j).Aff.RangeDiag);
elseif inp == 23
  M(k,i1:i2,m)=full(Prop{k,MatNo}(j).Asff.RangeDiag);
elseif inp == 24
  M(k,i1:i2,m)=full(Prop{k,MatNo}(j).Aff.RangeOffDiag);
elseif inp == 25
  M(k,i1:i2,m)=full(Prop{k,MatNo}(j).Asff.RangeOffDiag);
elseif inp == 26
  M(k,i1:i2,m)=full(Prop{k,MatNo}(j).Afc.RangeOffDiag);
elseif inp == 27
  M(k,i1:i2,m)=full(Prop{k,MatNo}(j).Asfc.RangeOffDiag);
elseif inp == 28
  M(k,i1:i2,m)=full(Prop{k,MatNo}(j).IterationTime);
elseif inp == 29
  M(k,i1:i2,m)=full(Prop{k,MatNo}(j).TimePerWayneOrderPerN);
elseif inp == 31
  M(k,i1:i2,m)=full(Prop{k,MatNo}(j).Aff.RangeDelta);
elseif inp == 32
  M(k,i1:i2,m)=full(Prop{k,MatNo}(j).Asff.RangeDelta);
 end;
end;
end;


for k=1:no
% number of levels for option k
 nl=size(Prop{k,MatNo},2);
 i1=0; i2=0;
 for j=1:nl;
  i1=i2+1; i2=i1+2;
if inp== 30;
  M(k,i1:i2,m)=full(Prop{k,MatNo}(j).WayneConvergence);
 end;
end;
end;
% end of numbers of matrices
end;

for m=1:TMatNo
  clear Mtable addtoheader headero;
  Mtable=M(:,:,m);
 disp('m; size of M(:,:,m)');disp(m);disp(size(Mtable));
  headero=char(Names(m,1));

%% Feed it to mTable2Tex
if inp <= 11 | inp==50 | inp==35 | inp==36;
  tit=[1:size(Mtable,2)];
elseif (inp >= 12 & inp < 32 & inp ~= 30)
  tit=[1:2*size(Mtable,2)]; 
elseif inp == 30
  tit=[1:3*size(Mtable,2)];
end;  
%
ctit=num2cell(tit);
disp('Hit any Key to Continue')
pause

iabbrev=0;
iabbrev=input(' Enter 1 to abbreviate titles: 0 Otherwise:');
fname=headero;
addtoheader =' ';
addtoheader = input('Enter description of table in quotes %s:')
header =[headero,'....', num2str(inp), '.... ', addtoheader]

mTable2Tex(Mtable,header,ctit,iabbrev,fname);
end;
 
