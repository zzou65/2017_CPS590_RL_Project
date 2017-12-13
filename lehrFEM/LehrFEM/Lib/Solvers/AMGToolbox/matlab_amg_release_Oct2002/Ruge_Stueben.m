%Ruge_Stueben	AMG interpolation matrix using Ruge-Stueben algorithm
% Modified nov 11 by jane

%usage : P=Ruge_Stueben(Aff,Afc,Astrongff,Astrongfc,opt)
%
%See also  AMGMakeInterpolation


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


function P=Ruge_Stueben(Aff,Afc,Astrongff,Astrongfc,opt)

  disp('   Ruge-Stueben')

  nc=size(Afc,2);
  nf=size(Aff,1);
  n=nc+nf;
  c=nf+[1:nc];
  f=1:nf;

  %% Do Ruge-Stueben stuff

  Af=[Aff,Afc];
  Astrongf=[Astrongff,Astrongfc];

  % Compensate for the weak elements
  if     strcmp( opt.weak, 'lump' )
    disp(['      Using lumped variant']);
    [Dff,Aff,Afc]=Lumping(Af,Astrongf,c,f,opt);
  elseif strcmp( opt.weak, 'scale' )
    disp(['      Using scaled variant']);
    [Dff,Aff,Afc]=Scaling(Af,Astrongf,c,f,opt);
  elseif strcmp( opt.weak, 'trash' )
    disp(['      Trashing weak elements']);
    % do nothing, just copy
    Dff=spdiag(diag(Aff));
    Aff=Astrongff;
    Afc=Astrongfc;
  else
    error(['      Unknown weak method ',opt.weak]);
  end

  % Make a guess for the unkown F-point errors (equation (6.4) stuff)
  [Dff,Afc]=Do64(Dff,Aff,Afc,opt);

  % Do equation 6.2
  if length(find(diag(Dff)==0))>0
    error(['Zeroes on diagonal of Dff, positions ',...
          num2str(find(diag(Dff)==0)')])
  end
  P=Dff\Afc;

return


%%%%%%%%%%

function [Dff,Aff,Afc]=Lumping(A,Astrong,c,f,opt)

  nc=length(c);
  nf=length(f);

  %% Do lumping (keep diagonal of Aff separate in Dff)

  Dff=spdiag(sum(A(f,:).')-sum(Astrong(f,:).'));
  Aff=Astrong(f,f);
  Afc=Astrong(f,c);

return


%%%%%%%%%%

function [Dff,Aff,Afc]=Scaling(A,Astrong,c,f,opt)

  nc=length(c);
  nf=length(f);
  n=nc+nf;

  %% Do scaling (keep diagonal of Aff separate in Dff)

  D=[spdiag(diag(A(f,f))),sparse(nf,nc)];

  Af=A(f,:)-D(f,:);
  Af(find(Af<0))=0;
  SumPos=full(sum(Af.'));

  Af=A(f,:)-D(f,:);
  Af(find(Af>0))=0;
  SumNeg=full(sum(Af.'));

  Afp=Astrong(f,:);
  Afp(find(Afp<0))=0;
  SumStrongPos=full(sum(Afp.'));

  Afn=Astrong(f,:);
  Afn(find(Afn>0))=0;
  SumStrongNeg=full(sum(Afn.'));

  StrongPos=find(SumStrongPos);
  StrongNeg=find(SumStrongNeg);

  ScalePos=zeros(nf,1);
  ScaleNeg=zeros(nf,1);
  ScalePos(StrongPos)=SumPos(StrongPos)./SumStrongPos(StrongPos);
  ScaleNeg(StrongNeg)=SumNeg(StrongNeg)./SumStrongNeg(StrongNeg);
  
  Afp=spdiag(ScalePos)*Afp;
  Afn=spdiag(ScaleNeg)*Afn;

  Aff=Afp(:,f)+Afn(:,f);
  Afc=Afp(:,c)+Afn(:,c);
  Dff=spdiag(diag(A(f,f)));

return


%%%%%%%%%%%%%

function [Dff,Afc_new]=Do64(Dff,Aff,Afc,opt)

  nc=size(Afc,2);
  nf=size(Aff,1);
  n=nc+nf;

  AffT=Aff.';
  clear Aff
  AfcT=Afc.';
  clear Afc
  AfcT_new=AfcT;

  if strcmp( opt.method64, 'original' )
    %% Do original equation 6.4 correction
    disp(['      Using original (6.4) correction (no abs)']);
    for i=1:nf
      if opt.nofill64
        Ci=find(AfcT(:,i));
      else
        Ci=[1:nc];
      end
      Fi=find(AffT(:,i));
      for j=Fi'
        % Get common C points 
        k=Ci(find(AfcT(Ci,j)));
        sumAfcjk=sum(AfcT(k,j));
        if sumAfcjk==0
          error(['Zero denominator in original 64, [i,j,k]=',...
                 num2str([i,j,k(:)'])])
        end
        AfcT_new(k,i)=AfcT_new(k,i)+AffT(j,i)*AfcT(k,j)/sumAfcjk;
      end
    end

  elseif strcmp( opt.method64, 'original abs' )
    %% Do original equation 6.4 correction with abs in sum
    disp(['      Using original (6.4) correction with abs in denominator']);
    for i=1:nf
      if opt.nofill64
        Ci=find(AfcT(:,i));
      else
        Ci=[1:nc];
      end
      Fi=find(AffT(:,i));
      for j=Fi'
        % Get common C points 
        k=Ci(find(AfcT(Ci,j)));
        AfcT_new(k,i)=AfcT_new(k,i)+AffT(j,i)*AfcT(k,j)/sum(-abs(AfcT(k,j)));
      end
    end

  elseif strcmp( opt.method64, 'original abs abs' )
    %% Do original equation 6.4 correction with abs in sum and numerator
    disp(['      Using original (6.4) correction with abs in denominator and numerator']);
    for i=1:nf
      if opt.nofill64
        Ci=find(AfcT(:,i));
      else
        Ci=[1:nc];
      end
      Fi=find(AffT(:,i));
      for j=Fi'
        % Get common C points 
        k=Ci(find(AfcT(Ci,j)));
        AfcT_new(k,i)=AfcT_new(k,i)+AffT(j,i)*abs(AfcT(k,j))/sum(abs(AfcT(k,j)));
      end
    end

  elseif strncmp( opt.method64, 'new', 3 )
    %% Do new equation 6.4 correction
    type=str2num(opt.method64(4:end));
    if length(type)<1 | type<1 | type>2
      error(['      Unknown 64 method ',opt.method64]);
    end
    disp(['      Using new (6.4) correction type ',...
          num2str(type)]);

    for i=1:nf
      if opt.nofill64
        Ci=find(AfcT(:,i));
      else
        Ci=[1:nc];
      end
      Fi=find(AffT(:,i));
      for j=Fi'
% modification below to maintain row sums
        fapos=1.0;
        faneg=1.0;
%%%%
	% Get common C points 
	k=Ci(find(AfcT(Ci,j)));
	Cj=find(AfcT(:,j));
	kpos = k(find(AfcT(k,j)>0));
	kneg = k(find(AfcT(k,j)<0));
% modification below to maintain row sums
        if type == 1
%          if (length(kpos)== 0 & length(kneg)>0), do nothing
%          if (length(kneg)== 0 & length(kpos)>0), do nothing
           if (length(kpos)>0 & length(kneg)> 0); fapos=0.5; faneg=0.5;end;
           if (length(kpos)==0 & length(kneg)==0); 
             error([' Fi to Fj Connection but No Common Coarse Point ',opt.method64]);
           end;
        end;
%%%
	if length(kpos)>0
	  Cjpos=Cj(find(AfcT(Cj,j)>0));
	  SumPosCj=sum(AfcT(Cjpos,j));
	  SumPosk =sum(AfcT(kpos, j));
	  ScalePos=SumPosCj/SumPosk;
	  if type == 1
% modification below to maintain row sums
	    AfcT_new(kpos,i)=AfcT_new(kpos,i)+fapos*AffT(j,i)*AfcT(kpos,j)/SumPosk;
	  else
	    AfcT_new(kpos,i)=AfcT_new(kpos,i)-AffT(j,i)*ScalePos*...
                                            AfcT(kpos,j)/Dff(j,j);
	  end
	end
	kneg = k(find(AfcT(k,j)<0));
	if length(kneg)>0
	  Cjneg=Cj(find(AfcT(Cj,j)<0));
	  SumNegCj=sum(AfcT(Cjneg,j));
	  SumNegk =sum(AfcT(kneg, j));
	  ScaleNeg=SumNegCj/SumNegk;
	  if type == 1
% modification below to maintain row sums
	    AfcT_new(kneg,i)=AfcT_new(kneg,i)+faneg*AffT(j,i)*AfcT(kneg,j)/SumNegk;
	  else
	    AfcT_new(kneg,i)=AfcT_new(kneg,i)-AffT(j,i)*ScaleNeg*...
                                            AfcT(kneg,j)/Dff(j,j);
	  end
	end
      end
    end

  else
    error(['      Unknown 64 method ',opt.method64]);
  end

  Afc_new=AfcT_new.';

return
