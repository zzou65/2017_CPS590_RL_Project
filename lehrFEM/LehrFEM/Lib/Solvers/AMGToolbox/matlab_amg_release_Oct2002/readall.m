%readall		Read matrix output from lamg and store in AMG data structure
%
%usage : L=readall
%
%Reads all data (*.rua and *.vec) in the current directory and 
%if available the grid data two directories down/up (../../*.grid)



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



function L=readall()

files=dir('*.rua');
for i=1:length(files)
  name=files(i).name;
  varname=stripext(name);
  [varname,level]=stripnr(varname);
  if level>0
    eval(['L{',num2str(level),'}.',varname,'=loadhb(name).'';']);
  end
end

files=dir('*.vec');
for i=1:length(files)
  name=files(i).name;
  varname=stripext(name);
  [varname,level]=stripnr(varname);
  if level>0
    [varname,nr]=stripnr(varname);
    if nr>0
      varname=[varname,'(:,',num2str(nr),')'];
    end
    if length(varname)>=4 & varname(1:4)=='jac_'
      varname=['M',varname(5:end)];
      eval(['L{',num2str(level),'}.',varname,'=spdiag(load(name));']);
    else
      eval(['L{',num2str(level),'}.',varname,'=load(name);']);
    end
    if length(varname)>=7 & varname(1:7)=='GSperm_'
      perm=eval(['L{',num2str(level),'}.',varname]);
      varname=['M',varname(8:end)];
      M=inv(tril(L{level}.A(perm,perm),0));
      [dummy,perm_inv]=sort(perm);
      M=M(perm_inv,perm_inv);
      eval(['L{',num2str(level),'}.',varname,'=M;']);
    end
  end
end


files=dir('../../*.grid');
if length(files)>1
  warning('multiple grid files :');
  files.name
end
if length(files)>0
  L{1}.mesh=load(['../../',files(1).name]);
  n=sqrt(size(L{1}.mesh,1));
  [t,n]=stripnr(files(1).name);
  t=t(1:end-1);
  [t,m]=stripnr(t);
  if n>0 & m>0
    disp(['Loading mesh size ',num2str(m),'x',num2str(n)]);
    [L{1}.grid,L{1}.meshx,L{1}.meshy]=mesh2coord(L{1}.mesh,m,n);
    L=propagategrid(L);
  else
    disp(['Could not get mesh size from file name ',files(1).name]);
  end
end

L{1}.c=L{1}.perm(1:L{1}.nc);
L{1}.f=L{1}.perm(L{1}.nc+1:end);
[dummy,L{1}.invperm]=sort(L{1}.perm);

function base=stripext(file)
  dot=find(file=='.');
  base=file(1:dot(1));

function [name,nr]=stripnr(name)
  name_=[name,' '];
  numbers=(name_>='0')&(name_<='9');
  start=find(numbers(1:end-1)==0 & numbers(2:end)==1)+1;
  stop =find(numbers(1:end-1)==1 & numbers(2:end)==0);
  if length(start)>0
    nr=str2num(name(start(end):stop(end)));
    name=name(1:start(end)-1);
    name_=find(name~='_');
    name=name(1:name_(end));
  else
    nr=-1;
  end
  
