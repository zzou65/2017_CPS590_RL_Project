% driver routine for several diffusion convection discretizations
% see below
  clear Mesh
  % Initialize constants
  JIG=2;
  d1=1;                      %impact of supg modification
  d2=1;
  a=10^-6;               %amount of diffusivity
  NREFS =5;
  % select test code
  V_Handle=...
            @(x,varargin)[ones(size(x,1),1) ones(size(x,1),1)];
        
  QuadRule = P7O6(); 
  Mesh.Coordinates =[-1 -1; 1 -1; 1 1; -1 1]; 
  Mesh.Elements = [1 2 4;2 3 4];
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = [-1 -1 -2 -2];   % -1 Inflow
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  
  err=zeros(NREFS,4);
  h=zeros(NREFS,1);
  Dofs=zeros(NREFS,1);
  
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
   
   %Laplace
   A = assemMat_LFE(NewMesh,@STIMA_Lapl_LFE);
   %Mass
   M_b=assemMat_LFE(NewMesh, @MASS_LFEquad, P1O2());
   M_v=assemMat_LFE(NewMesh, @MASS_LFEquad, P3O2());
   M_m=assemMat_LFE(NewMesh, @MASS_LFEquad, P3O3());
  
   %Convection term (upwinding)
   %diagonal Mass Matrix
   
   % upwinding Quadrature Rule= discrete differential forms
   I_1=assemMat_LFE(NewMesh, @STIMA_ContrGrad_Up, V_Handle);
   D_1b=M_b*I_1;
   D_1v=M_v*I_1;
   D_1m=M_m*I_1;
      
   % standard Galerkin D_0b=D0c
   D_0b=assemMat_LFE(NewMesh,@STIMA_Conv_LFE ,V_Handle,P1O2());
    
   % SUPG
   % D_supg1 = assemMat_LFE(NewMesh,@STIMA_Conv_LFE ,d.V_Handle,P1O2());
   D_supg2 = assemMat_LFE(NewMesh,@STIMA_SUPG_LFE,P1O2(), V_Handle,a,d1,d2);
   D_supg  = D_0b+D_supg2;
   
   % self-cooked
   M_QFELFE_b=assemMat_MASS_QFELFE(NewMesh,P1O2());
   M_QFELFE_v=assemMat_MASS_QFELFE(NewMesh,P3O2());
   M_QFELFE_m=assemMat_MASS_QFELFE(NewMesh,P3O3());
   M_QFELFE_e=assemMat_MASS_QFELFE(NewMesh,P7O4());
   %ContrOneQFE=assemMat_Contr1f_QFE(NewMesh,d.V_Handle);  % contraction of one forms
   %D_c=M_QFELFE*ContrOneQFE*TopGrad;
   I_2=assemMat_QFELFE(NewMesh,@STIMA_UPLFE_QFE,V_Handle);
   D_2b=M_QFELFE_b*I_2';
   D_2v=M_QFELFE_v*I_2';
   D_2m=M_QFELFE_m*I_2';
   D_2e=M_QFELFE_e*I_2';

   %Laplace + convection
  
   A_0b =(a*A+D_0b);
   A_1b =(a*A+D_1b);
   A_1v =(a*A+D_1v);
   A_1m =(a*A+D_1m);
   A_2b =(a*A+D_2b);
   A_2v =(a*A+D_2v);
   A_2m =(a*A+D_2m);
   A_2e =(a*A+D_2e);
   A_supg =(a*A+D_supg);
   
   % discrete inf-sup-condition
   
   FreeDofs = FreeDofs_LFE(NewMesh,[-1,-2]);
   B=assemMat_LFE_BdLayer(NewMesh,@STIMA_InfSup_LFE,QuadRule, V_Handle,1,0,0);

   E_0b = M_m(FreeDofs,FreeDofs)\A_0b(FreeDofs,FreeDofs); E_0b=A_0b(FreeDofs,FreeDofs)'*E_0b;
   
   E_1b = M_m(FreeDofs,FreeDofs)\A_1b(FreeDofs,FreeDofs); E_1b=A_1b(FreeDofs,FreeDofs)'*E_1b;
   E_1v = M_m(FreeDofs,FreeDofs)\A_1v(FreeDofs,FreeDofs); E_1v=A_1v(FreeDofs,FreeDofs)'*E_1v;
   E_1m = M_m(FreeDofs,FreeDofs)\A_1m(FreeDofs,FreeDofs); E_1m=A_1m(FreeDofs,FreeDofs)'*E_1m;
   
   E_2b = M_m(FreeDofs,FreeDofs)\A_2b(FreeDofs,FreeDofs); E_2b=A_2b(FreeDofs,FreeDofs)'*E_2b;
   E_2v = M_m(FreeDofs,FreeDofs)\A_2v(FreeDofs,FreeDofs); E_2v=A_2v(FreeDofs,FreeDofs)'*E_2v;
   E_2m = M_m(FreeDofs,FreeDofs)\A_2m(FreeDofs,FreeDofs); E_2m=A_2m(FreeDofs,FreeDofs)'*E_2m;
   E_2e = M_m(FreeDofs,FreeDofs)\A_2e(FreeDofs,FreeDofs); E_2e=A_2e(FreeDofs,FreeDofs)'*E_2e;
   
   E_supg = M_m(FreeDofs,FreeDofs)\A_supg(FreeDofs,FreeDofs); E_supg=A_supg(FreeDofs,FreeDofs)'*E_supg;
   
   err(i,1) = 1/eigs(inv(full(E_0b))*B(FreeDofs,FreeDofs),[],1,'LM');
   
   err(i,2) = 1/eigs(inv(full(E_1b))*B(FreeDofs,FreeDofs),[],1,'LM');
   err(i,3) = 1/eigs(inv(full(E_1v))*B(FreeDofs,FreeDofs),[],1,'LM');
   err(i,4) = 1/eigs(inv(full(E_1m))*B(FreeDofs,FreeDofs),[],1,'LM');
   
   err(i,5) = 1/eigs(inv(full(E_2b))*B(FreeDofs,FreeDofs),[],1,'LM');
   err(i,6) = 1/eigs(inv(full(E_2v))*B(FreeDofs,FreeDofs),[],1,'LM');
   err(i,7) = 1/eigs(inv(full(E_2m))*B(FreeDofs,FreeDofs),[],1,'LM');
   err(i,8) = 1/eigs(inv(full(E_2e))*B(FreeDofs,FreeDofs),[],1,'LM');
      
   err(i,9) = 1/eigs(inv(full(E_supg))*B(FreeDofs,FreeDofs),[],1,'LM');
  end; 
  tab=[1,2,3,4,9];
  
  fig = figure('Name','discr infsup 0 forms');
  plot(h,err(:,tab),'o-'); grid('on');
  set(gca,'XScale','log','YScale','log');
  Markers = '.ox+*sdv^<>ph';
  noMarkers = length(Markers);
  H1c = findobj(gca,'Type','line');
  linehandles = [H1c];
  for K = 1:length(linehandles)
     set(linehandles(K),'Marker',Markers(1+mod(K-1,noMarkers)));
  end
  xlabel(['{\bf h, diff = ',num2str(a),'}']);
  ylabel('{\bf infsup}');
  legend('(0,b)','(1,b)','(1,v)','(1,m)','(SUPG)','Location','Southeast')
  
  print -depsc './plots/infsup_LFEI1.eps'
  
  tab=[5,6,7,8,9];
  fig = figure('Name','discr infsup 0 forms');
  plot(h,err(:,tab),'o-'); grid('on');
  set(gca,'XScale','log','YScale','log');
  noMarkers = length(Markers);
  H1c = findobj(gca,'Type','line');
  linehandles = [H1c];
  for K = 1:length(linehandles)
     set(linehandles(K),'Marker',Markers(1+mod(K-1,noMarkers)));
  end
  xlabel(['{\bf h, diff = ',num2str(a),'}']);
  ylabel('{\bf infsup}');
  legend('(2,b)','(2,v)','(2,m)','(2,e)','(SUPG)','Location','Southeast')
  
  print -depsc './plots/infsup_LFEI2.eps'
  
clear all;