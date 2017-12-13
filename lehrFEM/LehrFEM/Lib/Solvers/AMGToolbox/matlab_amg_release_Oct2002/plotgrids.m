%plotgrids	Plot AMG grid points for the different levels
%
%usage : plotgrids(L[,mins[,maxs]])
%
%L : AMG data structure from AMGSetup. Must contain grid data.
%mins : Minimum marker size for gridpoint of the different levels,
%       default=8
%maxs : Maximum marker size, default=12
%
%See also  AMGPlotGrid


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



function plotgrids(L,mins,maxs)

if nargin<2
  mins=8;
end
if nargin<3
  maxs=12;
end

washold=ishold;

%colormap(gray);
colormap(jet);
c=colormap;
c=c(round(end/6):round(5*end/6),:);
levels=max(2,length(L));

for i=1:length(L)
%for i=levels:-1:1
  cnr=round((i-1)/(levels-1)*(length(c)-1)+1);
  size=(((i-1)*sqrt(maxs)+(levels-i)*sqrt(mins))/(levels-1))^2;
  plot(L{i}.grid(:,1),L{i}.grid(:,2),'d','markersize',size,...
       'MarkerEdgeColor','none','MarkerFaceColor',c(cnr,:));
  hold on;
end

if ~washold
  hold off
end
