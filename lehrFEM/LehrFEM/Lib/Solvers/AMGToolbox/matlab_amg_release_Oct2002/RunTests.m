%RunTests	Run tests with different matrices and options
%
%usage : RunTests([Options[,Names[,DoEig[,FileName]]]])
%
%This runs the AMGSetup and ComputeProperties for lots of option/matrix
%combinations and stores the results in an file
%
%Options  = Array of options structures, default = AMGDefaultOptions result
%Names    = cell array of matrix names, can contain wildcards, default='*.rua'
%DoEig    = if 1, will do the eigenvalue convergence factor in the
%           ComputeProperties, default=[] which means ComputeProperties will 
%           determine this using the size of the matrix 
%FileName = filename to store results in, default='Results_<date>_<time>.mat'
%
%Note: After each matrix and option combination the data is written to
%      a file labled "_incomplete", just in case something goes terribly
%      wrong. This file is removed after the final file is written.


% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% These programs were prepared by the Regents of the University of
% California at Los Alamos National Laboratory (the University) under
% contract W-7405-ENG-36 with the U.S. Department of Energy (DOE) and
% under Contract KC-07-01-01 with the Department of Energy(DOE),
% Office of Science, Mathematical, Information, and Computational Sciences,
% Applied Mathematical Sciences Program.  All rights in these programs
% are reserved by the DOE and the University.

% Permission is granted to the public to copy and use this software
% without charge, provided that this Notice and the statements of 
% authorship are reproduced on all copies. Neither the U.S. government 
% nor the University makes any warranty, express or implied, or assumes 
% any liability or responsibility for the use of this software.

% AMG code written by Menno Verbeek, in collaboration with Jane Cullum
%  and Wayne Joubert.
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


function RunTests(Options,Names,DoEig,filename)

if nargin<1 || isempty(Options)
  Options=AMGDefaultOptions;
end
if nargin<2 || isempty(Names)
  Names={'*.rua'};
end
if nargin<3 || isempty(DoEig)
  DoEig=[];
end
if nargin<4 || isempty(filename)
  filename=['Results ',datestr(now)];
end
filename=strrep(filename,' ','_');
filename=strrep(filename,':','_');
filename=strrep(filename,'.mat','');
tmpfilename=[filename,'_incomplete'];

if ~iscell(Names)
  Names={Names};
end
N={};
for i=1:length(Names)
  directory=dir(Names{i});
  N={N{:},directory.name};
end
Names=sort(N);
Names=Names(end:-1:1);
clear directory N

Names
filename

et.elapsed=[NaN,NaN];
et.nintervals=0;
ErrorTime.Total_Setup=et;
ErrorTime.Coarse_grid_selection=et;
ErrorTime.Pre_Smoother=et;
ErrorTime.Post_Smoother=et;
ErrorTime.Interpolation=et;
ErrorTime.Restriction=et;
ErrorTime.Coarse_grid_matrix=et;
ErrorTime.LU=et;
ErrorTime.Properties=et;

for i=1:length(Names)
  fprintf('***** %s (%d of %d) *****\n',Names{i},i,length(Names));
  %A=readhb(Names{i});
  A = hb_to_msm(Names{i});
  for j=1:length(Options)
    fprintf('*+*+* %s (%d of %d) *+*+*\n',Names{i},i,length(Names));
    fprintf('+++++ options %d of %d +++++\n',j,length(Options));
    lasterr('');
    try
      L=AMGSetup(A,Options(j));
      disp('Properties');
      starttimer('Properties');
      Prop{j,i}=ComputeProperties(L,DoEig);
      [t,f]=stoptimer('Properties');
      disp(['   Time = ',num2str(t),' MFlops = ',num2str(f/1e6)]);
      Times(j,i)=getalltimers;
      if isfield(Prop{j,i}(1),'Convergence') && ...
                      length(Prop{j,i}(1).Convergence)>0
        Conv(j,i)=Prop{j,i}(1).Convergence(1);
        disp(['   Convergence = ',num2str(Conv(j,i))]);
      else
        Conv(j,i)=Prop{j,i}(1).WayneConvergence(1);
        disp(['   WayneConvergence = ',num2str(Conv(j,i))]);
      end
      Opt(j)=L{1}.opt;
      %if ~comparestruct(mergestruct(AMGDefaultOptions,Options(j)),Opt(j))
      %  error('Options not consistent');
      %end
      Errors{j,i}='';
      clear L
    catch
      err=lasterr
      if strcmp( err, '' )
        return
      end
      Errors{j,i}=err;
      Times(j,i)=ErrorTime;
      Conv(j,i)=NaN;
      Prop{j,i}={};
      Opt(j)=mergestruct(AMGDefaultOptions,Options(j));
    end
    eval(['save ',tmpfilename,' Names Options Conv Times Opt Prop Errors']);
  end
end

eval(['save ',filename,' Names Options Conv Times Opt Prop Errors']);
disp('READY !!')
eval(['delete ',tmpfilename,'.mat'])
