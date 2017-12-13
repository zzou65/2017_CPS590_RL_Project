%showlocal	Attempt to plot local matrix connections in graph
%
%usage : showlocal(L,i,t)



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



function showlocal(L,i,t)

getstruct(L{1});

ip=find(perm==i);

if ip <= nc
  course=1;
  fine=0;
  disp('Course point');
else
  course=0;
  fine=1;
  disp('Fine point');
end

A2=A*A;
Ia=find(A(:,i));
d=find(Ia==i);
Ia(d)=Ia(1);
Ia(1)=i;

a=A(Ia,Ia);
as=Astrong(Ia,Ia);
s=S(Ia,Ia);
s(find(abs(s)<1e-10))=0;
a(find(abs(a)<1e-10))=0;
if nargin>2
  t=t(Ia);
end

figure;
showmat([a,zeros(size(a,1),1),as])
title('a ; as')

figure;
showmat(s)
title('s')

xy=[0,0];
dphi=2*pi/(ceil((length(Ia)-1)/2)*2+1);
for i=2:length(Ia)
  phi=(i-1+0.4*(rand-0.5))*dphi;
  xy(i,:)=[sin(phi),cos(phi)];
end

figure
hold on
for i=1:length(Ia)
  if find(perm==Ia(i))<=nc
    plot(xy(i,1),xy(i,2),'go','markersize',15)
  else
    plot(xy(i,1),xy(i,2),'ro','markersize',15)
  end
end
gplotval(as+diag(diag(a)),xy);
title('as+d')

figure
hold on
for i=1:length(Ia)
  if find(perm==Ia(i))<=nc
    plot(xy(i,1),xy(i,2),'go','markersize',15)
  else
    plot(xy(i,1),xy(i,2),'ro','markersize',15)
  end
end
gplotval(s,xy,.3);
title('s')

if nargin>2
figure
hold on
for i=1:length(Ia)
  if find(perm==Ia(i))<=nc
    plot(xy(i,1),xy(i,2),'go','markersize',15)
  else
    plot(xy(i,1),xy(i,2),'ro','markersize',15)
  end
end
gplotval(diag(t)+a-diag(diag(a)),xy,.3);
title('t')
end

Prow=P(i,:);
if course
  Pcol=P(:,ip);
    
else
  Pcol=[];
end



function gplotval(a,xy,f)

  if nargin<3
    f=0.5;
  end

  gplot(a,xy);
  hold on;
  for i=1:length(a)
    for j=1:length(a)
      if a(i,j)~=0 | i==j
        text(f*xy(i,1)+(1-f)*xy(j,1),f*xy(i,2)+(1-f)*xy(j,2),num2str(a(i,j),2));
      end
    end
  end

return
