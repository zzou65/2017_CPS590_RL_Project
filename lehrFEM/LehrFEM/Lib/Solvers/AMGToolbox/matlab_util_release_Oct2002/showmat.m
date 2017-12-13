%showmat		Show matrix using pcolor plot
%
%function showmat(A,climits)
%
%A matrix
%climits = [cmin,cmax] gives extremal matrix values to use for color map
%


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



function showmat(B,climits)

%clf;
[n,m]=size(B);

if ( length(find(imag(B)~=0)) > 0 )
  disp 'Using abs(A)'
  A=[full(abs(B)),zeros(n,1);zeros(1,m+1)];
else
  A=[full(real(B)),zeros(n,1);zeros(1,m+1)];
end
if (nargin<2)
  climits=[];
end

A(find(A==0))=NaN;

if ( size(climits,2) == 2 )
  h=climits(2);
  l=climits(1);
else
  h=max(max(A(1:end-1,1:end-1)));
  l=min(min(A(1:end-1,1:end-1)));
end
if l==h
   l=l-1;
   h=h+1;
end

climits=[l,h];
map=jet(128);
map=map(1:round(0.95*end),:);    % remove very dark red from map


%mesh(A);
%surf(A);
pcolor(A);

%axis('square');
axis('ij');
axis('image');
shading flat;
caxis(climits);
colormap(map);
colorbar

zoom on;

%colormap(cool)
%colormap(hot)
%colormap(gray)
%colormap(hsv)
%colormap(pink)
%colormap(copper)

hold off
%hold 
%hold
