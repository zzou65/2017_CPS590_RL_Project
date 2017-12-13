% driver routine for several diffusion convection discretizations
% see below
   clear Mesh
  % Initialize constants
  JIG=2;
  d1=1;                      %impact of supg modification
  d2=1;
  as=10.^-[4:10];               %amount of diffusivity
  NREFS =4;
  % select test code
  V_Handle=...
            @(x,varargin) [ones(size(x,1),1) ones(size(x,1),1)];
  %V_Handle=...
  %          @(x,varargin) [5+(x(:,1)+x(:,2)) 5+(x(:,1)+x(:,2))];
     
  QuadRule = P7O6(); 
   Mesh.Coordinates =[-1 -1; 1 -1; 1 1; -1 1]; 
  Mesh.Elements = [1 2 4;2 3 4];
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = [-1 -1 -2 -2];   % -1 Inflow
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  
  %err=zeros(NREFS,4);
  %h=zeros(NREFS,1);
  %Dofs=zeros(NREFS,1);
  
  for i = 1:NREFS
   
   %refine Mesh   
   Mesh = refine_REG(Mesh);
   
   % Mesh preprocessing
    
        switch(JIG)
        case 1
          NewMesh = Mesh;      
        case 2
          Loc = get_BdEdges(Mesh);
          Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
          FixedPos = zeros(size(Mesh.Coordinates,1),1);
          FixedPos(Loc) = 1;
          NewMesh = jiggle(Mesh,FixedPos);   
        case 3
          Loc = get_BdEdges(Mesh);
          Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
          FixedPos = zeros(size(Mesh.Coordinates,1),1);
          FixedPos(Loc) = 1;
          NewMesh = smooth(Mesh,FixedPos);
        end
   h(i)=get_MeshWidth(NewMesh);
   
  end
   
   %Laplace
   A = assemMat_QFE(NewMesh,@STIMA_Lapl_QFE);
   
   % linear interpolation
%    M_QFELFE_b=assemMat_MASS_QFELFE(NewMesh,P1O2());
%    M_QFELFE_v=assemMat_MASS_QFELFE(NewMesh,P3O2());
   M_QFELFE_m=assemMat_MASS_QFELFE(NewMesh,P3O3());
   M_QFELFE_e=assemMat_MASS_QFELFE(NewMesh,P4O3());
   
   I_1=assemMat_QFELFE(NewMesh, @STIMA_UPQFE_LFE, V_Handle);
   
%    D_1b=M_QFELFE_b'*I_1;
%    D_1v=M_QFELFE_v'*I_1;
   D_1m=M_QFELFE_m'*I_1;
   D_1e=M_QFELFE_e'*I_1;

   % standard Galerkin D_0b=D0c
   D_0b=assemMat_QFE(NewMesh,@STIMA_Conv_QFE ,V_Handle,P3O3());
     
   % quadratic interpolation
%    M_QFE_b=assemMat_QFE(NewMesh,@MASS_QFEquad,P1O2());
%    M_QFE_v=assemMat_QFE(NewMesh,@MASS_QFEquad,P3O2());
   M_QFE_vmb=assemMat_QFE(NewMesh,@MASS_QFEquad,P7O4());
   M_QFE_e=assemMat_QFE(NewMesh,@MASS_QFEquad,P7O6());
   M_QFE=assemMat_QFE(NewMesh,@MASS_QFEquad,P7O6());
   
   I_2=assemMat_QFE(NewMesh,@STIMA_UPQFE_QFE,V_Handle);
   
%    D_2b=M_QFE_b*I_2;
%    D_2v=M_QFE_v*I_2;
    D_2vmb=M_QFE_vmb*I_2;
    D_2e=M_QFE_e*I_2;
 
   [qr_bnd qr_inn]=split_QuadRule(P3O3());
   if (~isempty(qr_inn.x))
        B1=assemMat_QFE(NewMesh,@STIMA_Conv_QFE,V_Handle,qr_inn);
        if (~isempty(qr_bnd.x))
            B1=B1+assemMat_UPQFE2(NewMesh,@STIMA_UPQFE2,V_Handle,qr_bnd);
        end
   else
        if (~isempty(qr_bnd.x))
            B1=assemMat_UPQFE2(NewMesh,@STIMA_UPQFE2,V_Handle,qr_bnd);
        end
   end
   
   % weight2
   [qr_bnd qr_inn]=split_QuadRule(P4O3());
    if (~isempty(qr_inn.x))
        B2=assemMat_QFE(NewMesh,@STIMA_Conv_QFE,V_Handle,qr_inn);
        if (~isempty(qr_bnd.x))
            B2=B2+assemMat_UPQFE2(NewMesh,@STIMA_UPQFE2,V_Handle,qr_bnd);
        end
   else
        if (~isempty(qr_bnd.x))
            B2=assemMat_UPQFE2(NewMesh,@STIMA_UPQFE2,V_Handle,qr_bnd);
        end
   end
   
   % weight3
   [qr_bnd qr_inn]=split_QuadRule(P7O4());
    if (~isempty(qr_inn.x))
        B3=assemMat_QFE(NewMesh,@STIMA_Conv_QFE,V_Handle,qr_inn);
        if (~isempty(qr_bnd.x))
            B3=B3+assemMat_UPQFE2(NewMesh,@STIMA_UPQFE2,V_Handle,qr_bnd);
        end
    else
        if (~isempty(qr_bnd.x))
            B3=assemMat_UPQFE2(NewMesh,@STIMA_UPQFE2,V_Handle,qr_bnd);
        end
    end
   
   % discrete inf-sup-condition
   
   FreeDofs = FreeDofs_QFE(NewMesh,[-1,-2]);
   B=assemMat_QFE_BdLayer(NewMesh,@STIMA_InfSup_QFE,QuadRule, V_Handle,1,0,0);
   %B=assemMat_QFE(NewMesh,@STIMA_InfSup_QFE,QuadRule, V_Handle,1,0,0);
 
  %Laplace + convection
  for i=1:size(as,2)
        i 
        a=as(i);

       % SUPG
%        D_supg2 = assemMat_QFE(NewMesh,@STIMA_SUPG_QFE,P3O2(), V_Handle,a,d1,d2);
%        D_supg  = D_0b+D_supg2;
     %Laplace + convection
       A1=a*A+B1;
       A2=a*A+B2;
       A3=a*A+B3;
   
       A_0b =(a*A+D_0b);
%        A_1b =(a*A+D_1b);
%        A_1v =(a*A+D_1v);
       A_1m =(a*A+D_1m);
       A_1e =(a*A+D_1e);   
%        A_2b =(a*A+D_2b);
%        A_2v =(a*A+D_2v);
       A_2vmb =(a*A+D_2vmb);
       A_2e =(a*A+D_2e);
%       A_supg =(a*A+D_supg);
       

       % discrete inf-sup-condition
   
   %FreeDofs = FreeDofs_QFE(NewMesh,[-1,-2]);
   %B=assemMat_QFE_BdLayer(NewMesh,@STIMA_InfSup_QFE,QuadRule, V_Handle,1,0,0);
   %B=assemMat_QFE(NewMesh,@STIMA_InfSup_QFE,QuadRule, V_Handle,1,0,0);

   E_0b = M_QFE(FreeDofs,FreeDofs)\A_0b(FreeDofs,FreeDofs); E_0b=A_0b(FreeDofs,FreeDofs)'*E_0b;
   
%    E_1b = M_QFE(FreeDofs,FreeDofs)\A_1b(FreeDofs,FreeDofs); E_1b=A_1b(FreeDofs,FreeDofs)'*E_1b;
%    E_1v = M_QFE(FreeDofs,FreeDofs)\A_1v(FreeDofs,FreeDofs); E_1v=A_1v(FreeDofs,FreeDofs)'*E_1v;
   E_1m = M_QFE(FreeDofs,FreeDofs)\A_1m(FreeDofs,FreeDofs); E_1m=A_1m(FreeDofs,FreeDofs)'*E_1m;
   E_1e = M_QFE(FreeDofs,FreeDofs)\A_1e(FreeDofs,FreeDofs); E_1e=A_1e(FreeDofs,FreeDofs)'*E_1e;
   
%    E_2b = M_QFE(FreeDofs,FreeDofs)\A_2b(FreeDofs,FreeDofs); E_2b=A_2b(FreeDofs,FreeDofs)'*E_2b;
%    E_2v = M_QFE(FreeDofs,FreeDofs)\A_2v(FreeDofs,FreeDofs); E_2v=A_2v(FreeDofs,FreeDofs)'*E_2v;
   E_2vmb = M_QFE(FreeDofs,FreeDofs)\A_2vmb(FreeDofs,FreeDofs); E_2vmb=A_2vmb(FreeDofs,FreeDofs)'*E_2vmb;
   E_2e = M_QFE(FreeDofs,FreeDofs)\A_2e(FreeDofs,FreeDofs); E_2e=A_2e(FreeDofs,FreeDofs)'*E_2e;
  
   E1 = M_QFE(FreeDofs,FreeDofs)\A1(FreeDofs,FreeDofs); E1=A1(FreeDofs,FreeDofs)'*E1;
   E2 = M_QFE(FreeDofs,FreeDofs)\A2(FreeDofs,FreeDofs); E2=A2(FreeDofs,FreeDofs)'*E2;
   E3 = M_QFE(FreeDofs,FreeDofs)\A3(FreeDofs,FreeDofs); E3=A3(FreeDofs,FreeDofs)'*E3;
   
 %  E_supg = M_QFE(FreeDofs,FreeDofs)\A_supg(FreeDofs,FreeDofs); E_supg=A_supg(FreeDofs,FreeDofs)'*E_supg;
   
 err(i,1) = 1/max(eig(inv(full(E_0b))*full(B(FreeDofs,FreeDofs))));
   
%    err(i,2) = 1/eigs(inv(full(E_1b))*B(FreeDofs,FreeDofs),[],1,'LM');
%    err(i,3) = 1/eigs(inv(full(E_1v))*B(FreeDofs,FreeDofs),[],1,'LM');
   err(i,2) = 1/max(eig(inv(full(E_1m))*full(B(FreeDofs,FreeDofs))));
   err(i,3) = 1/max(eig(inv(full(E_1e))*full(B(FreeDofs,FreeDofs))));
   
%    err(i,5) = 1/eigs(inv(full(E_2b))*B(FreeDofs,FreeDofs),[],1,'LM');
%    err(i,6) = 1/eigs(inv(full(E_2v))*B(FreeDofs,FreeDofs),[],1,'LM');
%   err(i,5) = 1/eigs(inv(full(E_2m))*B(FreeDofs,FreeDofs),[],1,'LM');
   err(i,4) = 1/max(eig(inv(full(E_2vmb))*full(B(FreeDofs,FreeDofs))));
   err(i,5) = 1/max(eig(inv(full(E_2e))*full(B(FreeDofs,FreeDofs))));

   err(i,6) = 1/max(eig(inv(full(E1))*full(B(FreeDofs,FreeDofs))));
   err(i,7) = 1/max(eig(inv(full(E2))*full(B(FreeDofs,FreeDofs))));
   err(i,8) = 1/max(eig(inv(full(E3))*full(B(FreeDofs,FreeDofs))));
%    err(i,9) = 1/eigs(inv(full(E_supg))*B(FreeDofs,FreeDofs),[],1,'LM');
  end; 
  H=h(NREFS);
  tab=1:8;
  
  fig = figure('Name','discr infsup 0 forms');
  plot(as,err(:,tab),'o-'); grid('on');
  set(gca,'XScale','log','YScale','log');
  Markers = '.ox+*sdv^<>ph';
  noMarkers = length(Markers);
  H1c = findobj(gca,'Type','line');
  linehandles = [H1c];
  for K = 1:length(linehandles)
     set(linehandles(K),'Marker',Markers(1+mod(K-1,noMarkers)));
  end
  xlabel(['{\bf diff, h = ',num2str(H),'}']);
  ylabel('{\bf infsup}');
  legend('stand','(1,m)','(1,vb)','(2,vmb)','(2,ex)','(2,m)','(1b,vb)','(2b,vmb)','Location','NorthWest');
  print -depsc './plots/infsup_QFE_diff.eps'
  
%   tab=[5,6,7,8,9];
%   fig = figure('Name','discr infsup 0 forms');
%   plot(as,err(:,tab)); grid('on');
%   set(gca,'XScale','log','YScale','log');
%   xlabel(['{\bf diff, h = ',num2str(H),'}']);
%   ylabel('{\bf infsup}');
%   legend('(2,b)','(2,v)','(2,m)','(2,e)','(SUPG)','Location','SouthEast')
%   
%   print -deps 'infsup0formsI2_eps_qfe.eps'
  
clear all;