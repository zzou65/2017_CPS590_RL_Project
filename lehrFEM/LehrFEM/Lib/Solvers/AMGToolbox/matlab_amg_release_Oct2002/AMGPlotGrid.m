%AMGPlotGird	Plot the mesh, coarse/fine split and matrix connectivity
%
%function AMGPlotGird(L[,lstart[,lend]])
%
%L      : AMG data structure from AMGSetup, must contain the grid data
%lstart : first level to plot grid for, default 1
%lend   : last level to plot grid for, defaults to the lowest level in L


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


function AMGPlotGird(L,lstart,lend)

if nargin<2
  lstart=1;
end
if nargin<3
  lend=length(L);
end

washold=ishold;

if length(L{1}.meshx)>0
  plotmesh(L{1}.meshx,L{1}.meshy,'m-');
%  plotmesh(L{1}.meshx,L{1}.meshy,'g-');
  hold on;
end

if lstart<lend
  %gplot(L{lstart}.A,L{lstart}.grid);
  hold on
  gplot(L{lstart}.Astrong,L{lstart}.grid,'k');
end

if lend-lstart+1 > 2
  plotgrids({L{lstart:lend}});
else
  plotgrids({L{lstart:lend}},8,10);
end

zoom on;

if ~washold
  hold off
end


