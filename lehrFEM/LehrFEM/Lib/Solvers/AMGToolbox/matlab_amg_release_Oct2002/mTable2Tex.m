%Table2Tex	Make latex and .ps table output
%
%usage : mTable2Tex(Table[,header,Titles[,AbreviateTitles[,FileName]]])
%
%This will generate a table of the data in tex and ps format and
%show the ps file on screen using ghostview.
%
%Table : Cell or numeric array to make table of
%Titles: Cell array of strings with column titles, default: no titles
%AbreviateTitles : Whether to abreviate column titles (0 means full
%                  titles, 1 meand abreviated titles), default=0
%FileName : Filename base to use for the .tex and .ps file, default='table'
%
%See also  OptionsTable


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



function mTable2Tex(T,header,tit,abrvtit,name)

if isnumeric(T)
  T=num2cell(T);
end
if nargin<5
  name='table';
end
if nargin<4
  abrvtit=0;
end
if nargin<3 | length(tit)==0
  for j=1:size(T,2)
    tit{j}=['col(',num2str(j),')'];
  end
end

fid=fopen([name,'.tex'],'w');
%fid=1;

fprintf(fid,'\\documentclass{article}\n');
%fprintf(fid,'\\usepackage{lscape}\n');
%fprintf(fid,'\\setlength{\textheight}{23.3cm}\n');
fprintf(fid,'\\setlength{\\topmargin}{.01cm}\n');
fprintf(fid,'\\begin{document}\n');
%fprintf(fid,'\\landscape\n');
%fprintf(fid,'\\pagestyle{plain}\n');
%fprintf(fid,'\\small\n');
%fprintf(fid,'\\tiny\n');
fprintf(fid,'\\scriptsize\n');
%fprintf(fid,'\\footnotesize\n');
fprintf(fid,'\\noindent\n');

  fprintf(fid,'\\pagestyle{myheadings} \\markright{$');fprintf(fid,header);fprintf(fid,'$}');
  fprintf(fid,'\n');
%  fprintf(fid,'\\newpage');


%% Mark empty rows

for i=1:size(T,1)
  dorow(i)=0;
  for j=1:size(T,2)
    if length(T{i,j})>0 & ~strcmp(makestring(T{i,j}),'NaN')
      dorow(i)=1;
      break
    end
  end
end
Rows=find(dorow);
Rows=Rows(:)';

%% Remove constant columns

for j=1:size(T,2)
  keep(j)=0;
  for i=Rows
    if norm(size(T{Rows(1),j})~=size(T{i,j}),1) | ...
       ( length(T{Rows(1),j})>0 & norm(T{Rows(1),j}~=T{i,j},1) )
      keep(j)=1;
      break
    end
  end
end

%% Print constant columns

for j=find(keep==0)
  fprintf(fid,'%s = %s \\\\ \n',makestring(tit{j}),makestring(T{Rows(1),j}));
end

J=find(keep==1);

%% Make table

if abrvtit
  for j=1:length(J)
    coltit{J(j)}=j;
  end
else
  coltit=tit;
end

if abrvtit
  for j=J
    fprintf(fid,'col(%s) = %s \\\\ \n',makestring(coltit{j}),...
                makestring(tit{j}));
  end
end

%fprintf(fid,'\\begin{centeringn');
fprintf(fid,'\\hspace*{-4cm}\n');
fprintf(fid,'\\begin{tabular}{l|');
for j=J
  fprintf(fid,'l');
end
fprintf(fid,'}\n');

  fprintf(fid,' & ');
  for j=J
    fprintf(fid,' %s ',makestring(coltit{j}));
    if j<J(end)
      fprintf(fid,'&');
    end
  end
  fprintf(fid,'\\\\ \n');
  fprintf(fid,'\\hline \n');

for i=Rows
  fprintf(fid,' %d & ',i);
  for j=J
    fprintf(fid,' %s ',makestring(T{i,j}));
    if j<J(end)
      fprintf(fid,'&');
    end
  end
  fprintf(fid,'\\\\ \n');
end

fprintf(fid,'\\end{tabular}\n');
%fprintf(fid,'\\end{centering}');
fprintf(fid,'\\end{document}\n');

%!latex table.tex ; xdvi table &
%!latex table.tex ; dvips table -o ; ghostview table.ps &
eval(['!latex ',name,'.tex ; dvips ',name,' -o ; ',...
      'rm ',name,'.dvi ',name,'.log ',name,'.aux ; ',...
      'ghostview ',name,'.ps &']);
return
%%%%%%

function a=makestring(a)

  if ~ischar(a)
    a=num2str(a);
  else
    a=strrep(a,'\','$\backslash$');
    a=strrep(a,'_','\_');
    if length(a)>4 & a(end-3:end)=='.rua'
      a=a(1:end-4);
    end
    if 0
      if length(a)>=8 & a(1:8)=='augustus'
        a=['Aug',a(9:end)];
      end
      if length(a)>=7 & strcmpi( a(1:7), 'Poisson' )
        a=['Pois',a(8:end)];
      end
    end
  end

return
