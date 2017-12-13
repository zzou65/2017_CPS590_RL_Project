%AMGPlotVec	Plot a vector on the grid
%
%usage : AMGPlotVec(L,v)


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


function AMGPlotVec(L,v)

m=size(L{1}.meshx,1)-1;
n=size(L{1}.meshx,2)-1;
grid=L{1}.grid;

subplot(2,2,1)
draw(m,n,0,grid,v,'surfc');

subplot(2,2,2)
draw(m+1,n,m*n,grid,v,'surfc');

subplot(2,2,3)
draw(m+1,n,m*n+(m+1)*n,grid,v,'surfc');

subplot(2,2,4)
if 1
  cla;
  draw(m,n,0,grid,v,'mymesh');
  hold on;
  draw(m+1,n,m*n,grid,v,'mymesh');
  draw(m+1,n,m*n+(m+1)*n,grid,v,'mymesh');
  hold off;
else
  plotmesh(L{1}.meshx,L{1}.meshy,'m--');
  hold on;
  plotgrids(L,6,6);
  hold off;
end


function draw(m,n,offset,grid,v,type,varargin)

Gx=zeros(m,n);
Gy=zeros(m,n);
V=zeros(m,n);
Gx(:)=grid(offset+[1:m*n],1);
Gy(:)=grid(offset+[1:m*n],2);
V(:)=v(offset+[1:m*n]);

feval(type,Gx,Gy,V,varargin{:});


function mymesh(x,y,z,varargin)

surface(x,y,z,'FaceColor','none','EdgeColor','flat','Facelighting',...
        'none','EdgeLighting','Flat',varargin{:});
view(3);
grid on;
