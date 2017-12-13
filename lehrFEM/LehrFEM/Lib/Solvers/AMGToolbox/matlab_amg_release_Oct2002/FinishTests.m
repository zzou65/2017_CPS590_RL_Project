%FinishTests	Finish RunTests run or fill empty spots in RunTests data structure
%
%usage : FinishTests(R[,DoEig,[FileName]])
%
%This runs the AMGSetup and ComputeProperties for the open option/matrix 
%combinations in R and stores the new results in a file
%
%R        = Data structure of earlier RunTests run
%FileName = Filename to store results in, defaults to 
%           'Results_<date>_<time>.mat'


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


function RunTests(R,DoEig,filename)

if nargin<3 | length(DoEig)==0
  DoEig=[];
end
if nargin<3 | length(FileName)==0
  filename=['Results ',datestr(now)];
end
filename=strrep(filename,' ','_');

Options=R.Options
Names=R.Names
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
  fprintf('***** %s (%d of %d) *****',Names{i},i,length(Names));
  A=eval('readhb(Names{i})','[]');
  if length(A)>0
  for j=1:length(Options)
    fprintf('*+*+* %s (%d of %d) *+*+*\n',Names{i},i,length(Names));
    fprintf('+++++ options %d of %d +++++',j,length(Options));
    if length(R.Prop{j,i})>0
      disp('Already done');
    else
      lasterr('');
      try
        L=AMGSetup(A,Options(j));
        disp('Properties');
        starttimer('Properties');
        Prop{j,i}=ComputeProperties(L,DoEig);
        [t,f]=stoptimer('Properties');
        disp(['   Time = ',num2str(t),' MFlops = ',num2str(f/1e6)]);
        Times(j,i)=getalltimers;
        if isfield(Prop{j,i}(1),'Convergence') & ...
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
      eval(['save tmp_',filename,' Names Options Conv Times Opt Prop Errors']);
    end
  end
  end
end

eval(['save ',filename,' Names Options Conv Times Opt Prop Errors']);

