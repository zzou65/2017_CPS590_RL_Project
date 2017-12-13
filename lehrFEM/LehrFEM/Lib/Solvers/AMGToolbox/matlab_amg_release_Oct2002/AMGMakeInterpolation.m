%AMGMakeInterpolation	Make the interpolation matrix P
%
%usage : P=AMGMakeInterpolation(L,opt,transpose)


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


function P=AMGMakeInterpolation(L,opt,transpose)

n=length(L.perm);

%% No coarse points exception

if L.nc==0
  disp('   Warning : No coarse points');
  P=sparse(n,L.nc);
  return
end

%% No fine points exception

if L.nc==n
  disp('   Warning : No fine points');
  P=speye(n);
  return
end

%% Some setup

if ~strcmp( opt.matrix, 'A' ) && ...
   ~(strcmp( opt.objective, 'eigenvecs' ) && strcmp( opt.matrix, 'Se' ))
  if L.Mpre_explicit
    M=L.Mpre;
  else
    M=inv(L.Mpre);
  end
else
  M=[];
end

if transpose
  perm=L.permt;
  nc=L.nct;
  M=M.';
  A=L.A.';
  Astrong=L.Atstrong;
else
  perm=L.perm;
  nc=L.nc;
  %M=M;
  A=L.A;
  Astrong=L.Astrong;
end

n=length(perm);
c=perm(1:nc);
f=perm(nc+1:end);
nf=length(f);
[dummy,invperm]=sort(perm);

%% How to approximate it

% X=Minimize(A,B,As,Bs,C,opt) should minimize || AX-B || or X=A\B
% As,Bs give the strong variants of A and B
% C gives the (possible) sparsity of X, but might not be used
% opt gives the options

if     strcmp( opt.method, 'exact' )
  disp('   Using exact inverse/minimisation/eigenvecs/singularvecs')
  Approx=inline('A\B','A','B','As','Bs','C','opt');
elseif strcmp( opt.method, 'ainv' )
  disp(['   Using ainv, tol=',num2str(opt.ainv_tol)]);
  Approx=inline('ainv_explicit(A,1,opt.ainv_tol)*B',...
                  'A','B','As','Bs','C','opt');
elseif strcmp( opt.method, 'ainv strong' )
  disp(['   Using ainv on strong matrices, tol=',num2str(opt.ainv_tol)]);
  Approx=inline('ainv_explicit(As+spdiag(diag(A)),1,opt.ainv_tol)*Bs',...
                  'A','B','As','Bs','C','opt');
elseif strcmp( opt.method, 'gsai' )
  disp('   Using gsai')
  Approx=inline('gsai(A,B,C)','A','B','As','Bs','C','opt');
elseif strcmp( opt.method, 'Ruge-Stueben' ) || ...
       strcmp( opt.method, 'R-S' )
  disp('   Using (adapted) Ruge-Stueben method')
  Approx=inline('Ruge_Stueben(A,B,As,Bs,opt)','A','B','As','Bs','C','opt');
elseif strncmp( opt.method, 'local', 5 )
  disp(['   Using ',opt.method,' approximation method'])
  % Handled later
  Approx=inline('error(''   Option mismatch'')','A','B','As','Bs','C','opt');
else
  error(['   Unknown method ',opt.method]);
end


%% What sparsity pattern for the answer


if     strcmp( opt.pattern, 'dense' )
  disp('   Using dense pattern');
  Pattern=ones(nf,nc);
elseif strcmp( opt.pattern, 'Astrongfc' )
  disp('   Using pattern of Astrongfc');
  Pattern=spones(Astrong(f,c));
elseif strcmp( opt.pattern, 'Afc' )
  disp('   Using pattern of Afc');
  Pattern=spones(A(f,c));
elseif strcmp( opt.pattern, '' )
  disp('   Using empty pattern');
  Pattern=sparse(nf,nc);
else
  error(['   Unknown sparsity pattern ',opt.pattern]);
end


%% What matrix to use

if     strcmp( opt.matrix, 'A' )
  %% Standard
elseif strcmp( opt.matrix, 'MA' )
  %% Use MA in stead of A
  %% Do not use in combination with Ruge-Steuben or use_ainv==2
  disp('   Using M*A in stead of A (empty Astrong)');
  A=M*A;
elseif strcmp( opt.matrix, 'Se' )
  %% Use S=1-MA in stead of A
  %% Do not use in combination wiht Ruge-Steuben or use_ainv==2
  disp('   Using Se=1-M*A in stead of A (empty Astrong)');
  if strcmp( opt.objective, 'eigenvecs' ) && ~strcmp(opt.method, 'exact')
    global L_AMGBlindData
    L_AMGBlindData{1}=L;
    L_AMGBlindData{1}.opt.pre.its=1;
    A='AMGBlindS';
  else
    A=speye(n)-M*A;
  end
elseif strcmp( opt.matrix, 'Sr' )
  %% Use S=1-AM in stead of A
  %% Do not use in combination wiht Ruge-Steuben or use_ainv==2
  disp('   Using Sr=1-A*M in stead of A (empty Astrong)');
  A=speye(n)-A*M;
elseif strcmp( opt.matrix, 'Set' )
  %% Use S=1-M'A' in stead of A
  %% Do not use in combination wiht Ruge-Steuben or use_ainv==2
  disp('   Using Set=1-M''*A'' in stead of A (empty Astrong)');
  A=speye(n)-M'*A';
elseif strcmp( opt.matrix, 'Srt' )
  %% Use S=1-A'M' in stead of A
  %% Do not use in combination wiht Ruge-Steuben or use_ainv==2
  disp('   Using Srt=1-A''*M'' in stead of A (empty Astrong)');
  A=speye(n)-A'*M';
else
  error(['   Unknown matrix option ',opt.matrix]);
end


%% What to approximate

disp(['   Approximating ',opt.objective]);

if     strcmp( opt.objective, '0' )
  P=sparse(nf,nc);
  P=[speye(nc);P];
  P=P(invperm,:);
elseif strcmp( opt.objective, '-Aff\Afc' )
  P=-Approx(A(f,f),A(f,c),Astrong(f,f),Astrong(f,c),Pattern,opt);
  P=[speye(nc);P];
  P=P(invperm,:);
elseif strcmp( opt.objective, '-A:f\A:c' )
  P=-Approx(A(:,f),A(:,c),Astrong(:,f),Astrong(:,c),Pattern,opt);
  P=[speye(nc);P];
  P=P(invperm,:);
elseif strcmp( opt.objective, 'Aff/Acf' )
  % Same as (Acf'\Aff')'
  P=Approx(A(c,f).',A(f,f).',Astrong(c,f).',Astrong(f,f).',...
                   Pattern',opt).';
  P=[speye(nc);P];
  P=P(invperm,:);
elseif strcmp( opt.objective, 'Left-singular-vecs' )
  if strcmp(opt.method, 'exact')
    [u,s,v]=svd(full(A));
    P(:,1:nc)=u(:,1:nc);
  elseif strcmp(opt.method, 'local')
    [u,s,v]=svd(full(A));
    P(f,1:nc)=LocalSmoothApproximation(A(f,c),Astrong(f,c),u(c,1),u(f,1));
    P=[speye(nc);P];
    P=P(invperm,:);
  else
    error('   Option mismach')
  end
elseif strcmp( opt.objective, 'Right-singular-vecs' )
  if strcmp(opt.method, 'exact')
    [u,s,v]=svd(full(A));
    P(:,1:nc)=v(:,1:nc);
  elseif strcmp(opt.method, 'local')
    [u,s,v]=svd(full(A));
    P=LocalSmoothApproximation(A(f,c),Astrong(f,c),v(c,1),v(f,1));
    P=[speye(nc);P];
    P=P(invperm,:);
  else
    error('   Option mismach')
  end
elseif strcmp( opt.objective, 'eigenvecs' )
  if opt.eig_use_distance
    disp('   Using distances as weights');
    Astrong=GetDistances(L.grid,Astrong);
    I=find(Astrong);
    Astrong(I)=1./Astrong(I);
  end
  if strcmp(opt.method, 'exact')
    [e,ev]=sorteig(A);
    P(:,1:nc)=ev(:,end-nc+1:end);
  elseif strcmp(opt.method, 'local jdqr')
    [ev,e]=jdqr(A,1);
    if size(ev,2)==0
      error('JDQR did not converge');
    end
    ev=real(ev);
    assignin('caller','smooth',ev);
    P=LocalSmoothApproximation(Astrong(f,c),ev(c),ev(f));
    P=[speye(nc);P];
    P=P(invperm,:);
  elseif strcmp(opt.method, 'local jdqr new')
    [ev,e]=jdqr(A,1);
    if size(ev,2)==0
      error('JDQR did not converge');
    end
    ev=real(ev);
    assignin('caller','smooth',ev);
    P=LocalSmoothApproximationNew(Astrong(f,c),ev(c),ev(f));
    P=[speye(nc);P];
    P=P(invperm,:);
  elseif strcmp(opt.method, 'local') || strcmp(opt.method, 'local power')
    rand('state',0);
    ev=rand(n,1);
    ev=PowerMethod(A,ev,8);
    assignin('caller','smooth',ev);
    P=LocalSmoothApproximation(Astrong(f,c),ev(c),ev(f));
    P=[speye(nc);P];
    P=P(invperm,:);
  elseif strcmp(opt.method, 'local new') || strcmp(opt.method, 'local power new')
    rand('state',0);
    ev=rand(n,1);
    %ev=ones(n,1);
    ev=PowerMethod(A,ev,8);
    assignin('caller','smooth',ev);
    P=LocalSmoothApproximationNew(Astrong(f,c),ev(c),ev(f));
    P=[speye(nc);P];
    P=P(invperm,:);
  elseif strcmp(opt.method, 'local arnoldi new')
    rand('state',0);
    ev=rand(n,1);
    %ev=ones(n,1);
    ev=ArnoldiMethod(A,ev,8);
    assignin('caller','smooth',ev);
    P=LocalSmoothApproximationNew(Astrong(f,c),ev(c),ev(f));
    P=[speye(nc);P];
    P=P(invperm,:);
  elseif strcmp(opt.method, 'local 1')
    ev=ones(n,1);
    assignin('caller','smooth',ev);
    P=LocalSmoothApproximation(Astrong(f,c),ev(c),ev(f));
    P=[speye(nc);P];
    P=P(invperm,:);
  elseif strcmp(opt.method, 'local 1 new')
    ev=ones(n,1);
    assignin('caller','smooth',ev);
    P=LocalSmoothApproximationNew(Astrong(f,c),ev(c),ev(f));
    P=[speye(nc);P];
    P=P(invperm,:);
  elseif strcmp(opt.method, 'locallocal')
    P=LocalLocalSmoothApproximation(A,Astrong(f,c));
    P=[speye(nc);P];
    P=P(invperm,:);
  else
    error(['   Option mismach, method = ',opt.method])
  end
else
  error(['   Unknown interpolation objective ',opt.objective]);
end

%% Sparsify

if opt.droptol>0
  disp(['   Sparsifying P, tol=',num2str(opt.droptol)]);
  P(abs(P)<opt.droptol)=0;
end

%% Transpose ?

if transpose
  P=P.';
end

return


%%%%%%%%%%%%


function P=LocalSmoothApproximation(Astrongfc,smoothc,smoothf)

  nc=size(Astrongfc,2);
  nf=size(Astrongfc,1);
  n=nc+nf;

  P=sparse(nf,nc);
  for i=1:nf
    Ci=find(Astrongfc(i,:));
    sumabs=sum(abs(Astrongfc(i,Ci)));
    for j=Ci(:)'
      P(i,j)=smoothf(i)/smoothc(j)*abs(Astrongfc(i,j))/sumabs;
      %P(i,j)=smoothf(i)/smoothc(j)*1/length(Ci);
    end
  end

return


%%%%%%%%%%%%


function P=LocalSmoothApproximation2(Astrongfc,smoothc,smoothf)

  nc=size(Astrongfc,2);
  nf=size(Astrongfc,1);
  n=nc+nf;

  P=sparse(nf,nc);
  for j=1:nc
    Fj=find(Astrongfc(:,j));
    if length(Fj)>0
      P(Fj,j)=smoothf(Fj)/smoothc(j);
    end
  end
  for i=1:nf
    Ci=find(P(i,:));
    sumabs=sum(abs(Astrongfc(i,Ci)));
    P(i,Ci)=P(i,Ci).*abs(Astrongfc(i,Ci))/sumabs;
    %P(i,Ci)=P(i,Ci)*1/length(Ci);
  end

return

%%%%%%%%%%%%


function P=LocalSmoothApproximationNew(Astrongfc,smoothc,smoothf)

  nc=size(Astrongfc,2);
  nf=size(Astrongfc,1);
  n=nc+nf;
  AstrongfcT=Astrongfc.';

  PT=sparse(nc,nf);
  for i=1:nf
    Ci=find(AstrongfcT(:,i));
    if length(Ci)>0
      % Simple linear interpolation
      absAstrongfcTCii=abs(AstrongfcT(Ci,i));
      sumabs=sum(absAstrongfcTCii);
      PTCii=absAstrongfcTCii/sumabs;
      % Error in simple linear interpolation
      errori=smoothf(i)-PTCii.'*smoothc(Ci);
      % Correction
      sumabs=sum( abs(smoothc(Ci)) .* absAstrongfcTCii );
      PT(Ci,i)=PTCii+errori*sign(smoothc(Ci)).*absAstrongfcTCii/sumabs;
    end
  end

  P=PT.';

return

%%%%%%%%%%%%


function P=LocalLocalSmoothApproximation(S,Astrongfc)

  nc=size(Astrongfc,2);
  nf=size(Astrongfc,1);
  n=nc+nf;

  P=sparse(nf,nc);
  for j=1:nc
    smooth=S*S(:,j);
    Fj=find(Astrongfc(:,j));
    P(Fj,j)=smooth(nc+Fj)/smooth(j);
  end
  for i=1:nf
    Ci=find(P(i,:));
    if length(Ci)>0
      sumabs=sum(abs(Astrongfc(i,Ci)));
      P(i,Ci)=P(i,Ci).*abs(Astrongfc(i,Ci))/sumabs;
      %P(i,Ci)=P(i,Ci)*1/length(Ci);
    end
  end

return

%%%%%%%%%%%%


function v=PowerMethod(A,v,nr)

  for i=1:nr
    if ischar(A)
      v=eval([A,'(v)']);
    else
      v=A*v;
    end
    v=v/norm(v);
  end

return

%%%%%%%%%%%%


function v=ArnoldiMethod(A,v,nr)

  % Get Arnoldi AV~VH
  % Note: I get the square H
  %[H,V]=arnoldi(A,v,nr,0);
  [H,V]=myarnoldi(A,v,nr);
  [e,ev]=sorteig(H);
  v=V*ev(:,end);
  v=real(v);

return
