
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




function [time,fl,nintervals]=readtimer(name)

global internal_timer_data

if ~isfield(internal_timer_data,name)
  empty.running=0;
  empty.start=[0,0];
  empty.elapsed=[0,0];
  empty.nintervals=0;
  internal_timer_data=setfield(internal_timer_data,name,empty);
end

t=getfield(internal_timer_data,name);

if t.running
  warning('Timer is still running')
  time=t.elapsed(1)+cputime-t.start(1);
  flops=t.elapsed(2)+flops-t.start(2);
  nintervals=t.nintervals+0.5;
else
  time=t.elapsed(1);
  fl=t.elapsed(2);
  nintervals=t.nintervals;
end

