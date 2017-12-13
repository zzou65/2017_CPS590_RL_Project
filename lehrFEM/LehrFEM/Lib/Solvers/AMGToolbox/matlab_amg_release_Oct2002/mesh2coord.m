%mesh2coord	Convert mesh data to grid variable coordinates
%
%usage : [xy,Mx,My]=mesh2coord(mesh,m,n)
%
%See also  AMGSetup readall


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


function [xy,Mx,My]=mesh2coord(mesh,m,n)

if size(mesh,1)~=(m+1)*(n+1)
  error('Error, lenght mesh vectors doesn''t match');
end

dx=.1/m;
dy=.1/n;

Mx=ones(m+1,n+1);
My=Mx;
Mx(:)=mesh(:,1);
My(:)=mesh(:,2);

c=corner2center(Mx,My);
%c=c+[dx*ones(size(c,1),1),dy*ones(size(c,1),1)]
[f1,f2]=corner2face(Mx,My);

xy=[c;f1;f2];


return


%%%%%%%%%%

function centers=corner2center(x,y)

f=1;
x=(f*x(1:end-1,1:end-1)+x(2:end,1:end-1)+...
   x(1:end-1,2:end)+x(2:end,2:end))/(3+f);
y=(f*y(1:end-1,1:end-1)+y(2:end,1:end-1)+...
   y(1:end-1,2:end)+y(2:end,2:end))/(3+f);

centers=[x(:),y(:)];

%%%%%%%%%%

function [faces1,faces2]=corner2face(x,y)

x1=( x(1:end,1:end-1) + x(1:end,2:end) )/2;
y1=( y(1:end,1:end-1) + y(1:end,2:end) )/2;

x=x';
y=y';

x2=( x(1:end,1:end-1) + x(1:end,2:end) )/2;
y2=( y(1:end,1:end-1) + y(1:end,2:end) )/2;

faces1=[x1(:),y1(:)];
faces2=[x2(:),y2(:)];


