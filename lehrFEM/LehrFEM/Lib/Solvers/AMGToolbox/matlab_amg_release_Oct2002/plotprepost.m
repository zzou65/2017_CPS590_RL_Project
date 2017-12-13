%plotprepost	Do full eigen analysis of AMG problem
%
%This tool caluclates all kinds of diagnostics and plots them.
%
%usage : Lout=plotprepost(L[,level[,doall]])
%
%L     : AMG data structure form AMGSetup
%level : Level in the AMG structure to use, default=1
%doall : Recalculate everything if doall==1, default=0
%Lout  : a copy of L with all the diagnostic information added
%
%Note 1 : If an error or interrupt (^C) occurs, the calculated data will
%         be returned in Lout. Calculation can be resumed using
%         >> Lout=plotprepost(Lout,level)
%
%Note 2 : Figures can be replot without redoing the calculation using
%         >> plotprepost(Lout,level);


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


function L=plotprepost(L,l,doall)

if nargin<2
  l=1;
end
if nargin<3
  doall=0;
end

opt=L{1}.opt;

if norm(L{l}.Mpre-L{l}.Mpost,1)>0
  disp('Warning: Mpre~=Mpost, using Mpre only for both')
  disp('         This can be caused by resuming a plotprepost calculation')
end
if ~L{l}.Mpre_explicit
  L{l}.Mpre=inv(L{l}.Mpre);
  L{l}.Mpre_explicit=1;
end
if max(max(abs(L{l}.Ac-L{l}.R*L{l}.A*L{l}.P)))>1e-10
  disp(['Warning: Ac ~= R*A*P, max error is ',...
       num2str(max(max(abs(L{l}.Ac-L{l}.R*L{l}.A*L{l}.P))))]);
  if max(max(abs(L{l}.Ac-L{l}.R*L{l}.A*L{l}.P)))>1e-8
    disp('Please execute L{l}.Ac=L{l}.R*L{l}.A*L{l}.P');
    return;
  else
    disp('Doing Ac=R*A*P');
    L{l}.Ac=L{l}.R*L{l}.A*L{l}.P;
  end
end

try

%%%%%%%%%%%% Operator calculation

id=speye(size(L{l}.A));

%if ~isfield(L{l},'Acinv')
%  disp('Inverting Ac')
%  L{l}.Acinv=full(inv(L{l}.Ac));
%end

if ~isfield(L{l},'T') | doall
  disp('Calcluating G, S, K and T')
  L{l}.G=GetG(L,l);
  if L{l}.Mpre_explicit
    L{l}.S=(id-L{l}.Mpre*L{l}.A)^opt.pre.its;
  else
    L{l}.S=(id-L{l}.Mpre\L{l}.A)^opt.pre.its;
  end
  if nnz(L{l}.S)/prod(size(L{l}.S)) > 0.2
    L{l}.S=full(L{l}.S);
  end
  L{l}.K=id-L{l}.G*L{l}.A;
  L{l}.T=L{l}.S*L{l}.K*L{l}.S;
end


%%%%%%%%%%% eigenvalues/vectors calcualtion


if ~isfield(L{l},'et') | doall
  disp('Calculating max eigenvalue of T')
  conv=jdqr(L{l}.T,1)
  disp('Calculating eingvectors of T')
  [L{l}.et,L{l}.evt]=sorteig(L{l}.T,'real');
  if (max(abs(imag(L{l}.et))))<1e-8
    L{l}.et=real(L{l}.et);
  else
    disp('Warning: eigenvalues of T are complex')
  end
end
[dummy,Itmp]=sort(abs(L{l}.et));
tmp=L{l}.et(Itmp(end));
disp(['*** Maximum eigenvalue of S*K*S : | ',num2str(tmp),' | = ',...
     num2str(abs(tmp))]);
if ~isfield(L{l},'es') | doall
  disp('Calculating eingvectors of S')
  [L{l}.es,L{l}.evs]=sorteig(L{l}.S,'real');
  if (max(abs(imag(L{l}.es))))<1e-8
    L{l}.es=real(L{l}.es);
  else
    disp('Warning: eigenvalues of S are complex')
  end
end
if ~isfield(L{l},'ek') | doall
  disp('Calculating eingvectors of K')
  [L{l}.ek,L{l}.evk]=sorteig(L{l}.K,'real');
  if (max(abs(imag(L{l}.ek))))<1e-8
    L{l}.ek=real(L{l}.ek);
  else
    disp('Warning: eigenvalues of K are complex')
  end
end
if ~isfield(L{l},'ea') | doall
  disp('Calculating eingvectors of A')
  [L{l}.ea,L{l}.eva]=sorteig(L{l}.A,'real');
  if (max(abs(imag(L{l}.ea))))<1e-8
    L{l}.ea=real(L{l}.ea);
  else
    disp('Warning: eigenvalues of A are complex')
  end
end
if ~isfield(L{l},'Nr') | doall
  disp('Calculating null(R)')
  L{l}.Nr=null(full(L{l}.R));
end


%%%%%%%%%%%%%% norm reductions

% get rid of old fields
fn=fieldnames(L{l});
for i=1:length(fn)
  if length(fn{i})>2 & fn{i}(1)=='n'
    if size(getfield(L{l},fn{i}))==[size(L{l}.A,1),1]
      if ~sum(strcmp(fn{i},{'nsk','nske1','nsnr'}))
        disp(['Removing old field ',fn{i}]);
        L{l}=rmfield(L{l},fn{i});
      end
    end
  end
end


% effect of S*K*S on S*K*S eigenvectors

L{l}.nrt=plotnormred(L{l},L{l}.evt,L{l}.et,'S*K*S','nrt','abs(e)',doall);
L{l}.nrt=plotnormred(L{l},L{l}.evt,L{l}.et,'S*K*S','nrt','abs(nrs1)',doall);
L{l}.nrt=plotnormred(L{l},L{l}.evt,L{l}.et,'S*K*S','nrt','abs(nrk)',doall);
L{l}.nrt=plotnormred(L{l},L{l}.evt,L{l}.et,'S*K*S','nrt','abs(nrs2)',doall);


% effect of S*K*S on A eigenvectors

L{l}.nra=plotnormred(L{l},L{l}.eva,L{l}.ea,'A','nra','abs(e)',doall);
L{l}.nra=plotnormred(L{l},L{l}.eva,L{l}.ea,'A','nra','abs(nrt)',doall);


% norm reductions by S

L{l}.nrs=plotnormred(L{l},L{l}.evs,L{l}.es,'S','nrs','abs(e)',doall);
L{l}.nrs=plotnormred(L{l},L{l}.evs,L{l}.es,'S','nrs','abs(nrt)',doall);

if ~isfield(L{l},'nsk') | doall
  L{l}.nsk=norms(L{l}.S*L{l}.evk);
end
if ~isfield(L{l},'nske1') | doall
  Ike1=find(abs(L{l}.ek)>.5);
  L{l}.nske1=L{l}.nsk(Ike1);
end
if ~isfield(L{l},'nsnr') | doall
  L{l}.nsnr=norms(L{l}.S*L{l}.Nr);
end

getstruct(L{l});

figure
plot(abs(es),'g');
hold on
plot(sort(nrt.nrs1),'b')
plot(sort(nrt.nrs2),'k')
plot(sort(nsnr),'r')
plot(sort(nske1),'m')
title('norm reductions by S')
xlabel('Sorted by size')
legend('S eigenvalues','pre-smooth norm reduction in S*K*S',...
       'post-smooth norm reduction in S*K*S',...
       'norm reduction by S on null(R)',...
       'norm reduction by S on 1 eigenspace of K')

% norm reductions by K

figure
plot(abs(ek),'m');
hold on
plot(sort(nrt.nrk),'r')
plot(nsk,'k')
title('norm reductions by K')
xlabel('Sorted by K eigenvalues / size')
legend('K eigenvalues','norm reduction in S*K*S',...
       'norm reduction by S on eigenvectors of K')



%%%%%%%%% Conformation to eq. (6.2)

if 0
% get rid of old fields
fn=fieldnames(L{l});
for i=1:length(fn)
  if length(fn{i})>4 & fn{i}(1:3)=='c62'
    disp(['Removing old field ',fn{i}]);
    L{l}=rmfield(L{l},fn{i});
  end
end

L{l}.c62s=plot62conf(L{l},L{l}.evs,L{l}.es,'S','c62s',doall);
L{l}.c62a=plot62conf(L{l},L{l}.eva,L{l}.ea,'A','c62a',doall);
L{l}.c62k=plot62conf(L{l},L{l}.evk,L{l}.ek,'K','c62k',doall);
L{l}.c62p=plot62conf(L{l},L{l}.P,[],'P','c62p',doall);
end

%%%%%%%% Representation of ev's on P

disp('Representation of ev''s on P')
if ~isfield(L{l},'prs') | doall
  L{l}.prs=norms(L{l}.P*(full(L{l}.P)\L{l}.evs)-L{l}.evs);
end
if ~isfield(L{l},'rrs') | doall
  L{l}.rrs=norms(L{l}.R'*(full(L{l}.R')\L{l}.evs)-L{l}.evs);
end

figure
plot(L{l}.prs,'b');
hold on
plot(L{l}.rrs,'r');
plot(abs(L{l}.es),'g');
title('norms(evs-P*pinv(P)*evs)');
xlabel('Sorted by S eigenvalues')
legend('norms(evs-P*pinv(P)*evs)','norms(evs-R''*pinv(R'')*evs)','S eigenvalues')

if ~isfield(L{l},'pra') | doall
  L{l}.pra=norms(L{l}.P*(full(L{l}.P)\L{l}.eva)-L{l}.eva);
end
if ~isfield(L{l},'rra') | doall
  L{l}.rra=norms(L{l}.R'*(full(L{l}.R')\L{l}.eva)-L{l}.eva);
end

figure
plot(L{l}.pra,'b');
hold on
plot(L{l}.rra,'r');
plot(L{l}.ea/max(abs(L{l}.ea)),'g');
title('norms(eva-P*pinv(P)*eva)');
xlabel('Sorted by A eigenvalues')
legend('norms(eva-P*pinv(P)*eva)','norms(eva-R''*pinv(R'')*eva)',...
       ['scaled A eigenvalues ',num2str(max(abs(L{l}.ea)))])

%%%%%%%% Projections

if ~isfield(L{l},'Snorms') | doall
  disp('Calculating projections')
  c=L{l}.c;
  f=L{l}.f;
  nc=length(c);
  nf=length(f);
  n=nc+nf;
  cc=1:nc;
  ff=nc+1:n;
  PP(1:n,cc)=L{l}.P;
  PP(f,ff)=speye(nf);
  PP(c,ff)=0;
  RR(cc,1:n)=L{l}.R;
  RR(ff,f)=speye(nf);
  RR(ff,c)=0;
  L{l}.Arp=RR*L{l}.A*PP;
  L{l}.Spp=inv(PP)*L{l}.S*PP;
  L{l}.Kpp=inv(PP)*L{l}.K*PP;
  L{l}.Tpp=inv(PP)*L{l}.T*PP;

  L{l}.Sppnorms=[norm(full(L{l}.Spp(cc,cc))),norm(full(L{l}.Spp(cc,ff)));...
                 norm(full(L{l}.Spp(ff,cc))),norm(full(L{l}.Spp(ff,ff)))];
  L{l}.Kppnorms=[norm(full(L{l}.Kpp(cc,cc))),norm(full(L{l}.Kpp(cc,ff)));...
                 norm(full(L{l}.Kpp(ff,cc))),norm(full(L{l}.Kpp(ff,ff)))];
  L{l}.Tppnorms=[norm(full(L{l}.Tpp(cc,cc))),norm(full(L{l}.Tpp(cc,ff)));...
                 norm(full(L{l}.Tpp(ff,cc))),norm(full(L{l}.Tpp(ff,ff)))];
  L{l}.Snorms=[norm(full(L{l}.S(c,c))),norm(full(L{l}.S(c,f)));...
               norm(full(L{l}.S(f,c))),norm(full(L{l}.S(f,f)))];
end

figure
subplot(2,2,1)
showmat(L{l}.Snorms);
title('S block norms')
subplot(2,2,2)
showmat(L{l}.Sppnorms);
title('Spp block norms')
subplot(2,2,3)
showmat(L{l}.Kppnorms);
title('Kpp block norms')
subplot(2,2,4)
showmat(L{l}.Tppnorms);
title('Tpp block norms')

catch

  disp('An error or interrupt occurd, calculated data is returned');
  disp('To continue the calculation, do L=plotprepost(L);');

end


% end main

%%%%%%%%%%% Subroutines

function data=plotnormred(L,ev,e,name,varname,order,doall)
  version=1;

  if ~isfield(L,varname) | eval(['L.',varname,'.version'])~=version | doall
    disp(['Calculating the norm reductions for eigenvectors of ',name])
    data.version=version;

    v=ev;
    v=v*spdiag(1./norms(v));
    v=L.S*v;
    data.nrs1=norms(v);
    v=v*spdiag(1./data.nrs1);
    v=L.K*v;
    data.nrk=norms(v);
    v=v*spdiag(1./data.nrk);
    v=L.S*v;
    data.nrs2=norms(v);
    data.nrt=data.nrs1.*data.nrk.*data.nrs2;
  else
    data=eval(['L.',varname]);
  end

  getstruct(data);
  [dummy,I]=sort(eval(order));
  clear dummy;

  figure
  plot(nrt(I),'r');
  hold on
  plot(nrs1(I),'b');
  plot(nrk(I),'m');
  plot(nrs2(I),'k');
  plot(abs(e(I)),'g');
  plot(nrt(I),'r');
  title(['effect of S*K*S on ',name,' eigenvectors, max. nrt.=',num2str(max(nrt))])
  xlabel(['Sorted by ',order])
  legend('norm reduction by S*K*S','norm reduction by S pre',...
         'norm reduction by K','norm reduction by S post',...
         [name,' eigenvalues'])
  drawnow;

% end plotnormred


function data=plot62conf(L,ev,e,name,varname,doall)
  version=3;

  if ~isfield(L,varname) | eval(['L.',varname,'.version'])~=version | doall
    disp(['Calculating the conformation to eq. (6.2) for ',name])
    data.version=version;

    c=L.perm(1:L.nc);
    f=L.perm(L.nc+1:end);
    Ax=L.A*ev;
    Dx=spdiag(diag(L.A))*ev;
    data.normc=norms(Dx(c,:))./norms(Ax(c,:));
    data.normf=norms(Dx(f,:))./norms(Ax(f,:));
    Ax=full(Ax);
    Dx=full(Dx);
    nullen=find( (abs(Ax)<1e-10) & (abs(Dx)<1e-10) );
    Ax(nullen)=NaN;
    Dx(nullen)=NaN;
    tc=sort(abs(Dx(c,:)./Ax(c,:)));
    tf=sort(abs(Dx(f,:)./Ax(f,:)));
    for i=1:size(tc,2)
      uc=tc(find(~isnan(tc(:,i))),i);
      uf=tf(find(~isnan(tf(:,i))),i);
      %data.meanc(i)=mean(uc);
      %data.meanf(i)=mean(uf);
      if length(uc)>0
        data.medc(i)=uc(ceil(length(uc)/2));
      else
        data.medc(i)=NaN;
      end
      if length(uf)>0
        data.medf(i)=uf(ceil(length(uf)/2));
      else
        data.medf(i)=NaN;
      end
    end
  else
    data=eval(['L.',varname]);
  end

  figure
  plot(data.normf,'b');
  hold on
  plot(data.medf,'k');
  plot(data.normc,'r');
  plot(data.medc,'m');
  if length(e)>0
    plot(abs(e),'g');
  end
  title(['Conformation to eq. (6.2) of eigenvectors of ',...
        name,' (large is good)'])
  xlabel(['Sorted by ',name,' eigenvalue size']);
  legend('fine |Dx|/|Ax|','fine median of Dx./Ax',...
         'coarse |Dx|/|Ax|','coarse median of Dx./Ax',...
         [name,' eigenvalue']);
  a=axis;
  axis([a(1:2),-1,5])
  drawnow;

% end function plot62conf


function G=GetG(L,l)

G=L{l}.P*GetAcinv(L,l+1)*L{l}.R;


function t=GetAcinv(L,l)

if l<length(L)
  id=speye(size(L{l}.A,2));

  if L{l}.Mpre_explicit
    Mpre=L{l}.Mpre;
  else
    Mpre=inv(L{l}.Mpre);
  end
  if L{l}.Mpost_explicit
    Mpost=L{l}.Mpost;
  else
    Mpost=inv(L{l}.Mpost);
  end

  t=Mpre;
  for i=2:L{1}.opt.pre.its
    t=t+Mpre*(id-L{l}.A*t);
  end
  t=t+GetG(L,l)*(id-L{l}.A*t);
  for i=1:L{1}.opt.post.its
    t=t+Mpost*(id-L{l}.A*t);
  end
else
  t=inv(L{l}.U)*inv(L{l}.L);
end

