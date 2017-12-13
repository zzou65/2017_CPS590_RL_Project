%GetDistances	Get distance matrix from the grid (2D only)
%
%usage : D=GetDistances(L[,sparse[,level]]) or D=GetDistances(grid[,pattern])
%
%L     = AMGSetup data structure
%level = levet to use from L, default=1
%grid  = grid data from mesh2coord, must contain 2 columns giving
%        the x- and y-coordiate of the grid points
%sparse= if 0, D is dense, if 1, return D with sparsity pattern of A, default=1
%pattern = if set, use this sparsety pattern, default dense
%D     = symmetric dense matrix with D_ji the distance between grid points i and j


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


function D=GetDistances(L,sparse,l)

if nargin<3
  l=1;
end

if iscell(L) & isfield(L{l},'grid')
  grid=L{l}.grid;
  A=L{l}.A;
  if nargin<2
    sparse=1;
  end
else
  grid=L;
  if nargin<2
    sparse=0;
  else
    A=sparse;
    sparse=1;
  end
end

grid=grid(:,1)+i*grid(:,2);
n=length(grid);

if sparse
  [I,J]=find(A);
  for i=1:length(I)
    D(i,1)=abs(grid(I(i))-grid(J(i)));
  end
  D=spconvert([I,J,D]);
else
  D=abs(ones(n,1)*grid.'-grid*ones(1,n));
end
