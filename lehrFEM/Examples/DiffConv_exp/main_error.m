% driver routine for several diffusion convection discretizations
% see below
   clear Mesh
  % Initialize constants
  JIG=2;
  d1=1;                      %impact of supg modification
  d2=1;
  a=10^0;               %amount of diffusivity
  NREFS =7;
  % select test code
  d=getData(5);
     
  QuadRule = P7O6(); 
  
  Mesh.Coordinates =d.Coordinates; 
  Mesh.Elements = d.Elements;
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = d.boundtype;
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  
  err=zeros(NREFS,4);
  h=zeros(NREFS,1);
  l=zeros(NREFS,2);
  Dofs=zeros(NREFS,1);
  
  for i = 1:NREFS
   i
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
   I_1=assemMat_LFE(NewMesh, @STIMA_ContrGrad_Up, d.V_Handle);
   
%    LU = assemLoad_LFE(NewMesh,P7O6(),d.U_EX_Handle,a);
%    l(i,1)=L2Err_LFE(NewMesh,I_1*(M_m\LU),P7O6(),d.V_GRAD_U_EX_Handle,0,a)
   D_1b=M_b*I_1;
   D_1v=M_v*I_1;
   D_1m=M_m*I_1;
      
   % standard Galerkin D_0b=D0c
   D_0b=assemMat_LFE(NewMesh,@STIMA_Conv_LFE ,d.V_Handle,P1O2());
    
   % SUPG
   % D_supg1 = assemMat_LFE(NewMesh,@STIMA_Conv_LFE ,d.V_Handle,P1O2());
   D_supg2 = assemMat_LFE(NewMesh,@STIMA_SUPG_LFE,P1O2(), d.V_Handle,a,d1,d2);
   D_supg  = D_0b+D_supg2;
   
   % self-cooked
   M_QFELFE_b=assemMat_MASS_QFELFE(NewMesh,P1O2());
   M_QFELFE_v=assemMat_MASS_QFELFE(NewMesh,P3O2());
   M_QFELFE_m=assemMat_MASS_QFELFE(NewMesh,P3O3());
   M_QFELFE_e=assemMat_MASS_QFELFE(NewMesh,P7O4());
   %ContrOneQFE=assemMat_Contr1f_QFE(NewMesh,d.V_Handle);  % contraction of one forms
   %D_c=M_QFELFE*ContrOneQFE*TopGrad;
   I_2=assemMat_QFELFE(NewMesh,@STIMA_UPLFE_QFE,d.V_Handle);
   %LU = assemLoad_QFE(NewMesh,P7O6(),d.U_EX_Handle,a);
   %l(i,2)=L2Err_QFE(NewMesh,I_2'*(M_m\LU),P7O6(),d.V_GRAD_U_EX_Handle,0,a)
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
   % source term
  
   L = assemLoad_LFE(NewMesh,P7O6(),d.SOL_Handle,a);
   L_U = assemLoad_LFE(NewMesh,P7O6(),d.U_EX_Handle,a);
   L_supg=assemLoad_LFE_SUPG(NewMesh,P1O2(),d.V_Handle,d.SOL_Handle,a,d1,d2);
   
   % Direchlet boundary
  
   [U_0b,FreeDofs] = assemDir_LFE(NewMesh,-1,d.U_EX_Handle,a);
   U_1b=U_0b;
   U_1v=U_0b;
   U_1m=U_0b;
   U_2b=U_0b;
   U_2v=U_0b;
   U_2m=U_0b;
   U_2e=U_0b;
   U_supg=U_0b;
   
   L_0b = L - A_0b*U_0b;
   L_1b = L - A_1b*U_1b;
   L_1v = L - A_1v*U_1v;
   L_1m = L - A_1m*U_1m;
   L_2b = L - A_2b*U_2b;
   L_2v = L - A_2v*U_2v;
   L_2m = L - A_2m*U_2m;
   L_2e = L - A_2e*U_2e;

   L_supg = (L+L_supg) - A_supg*U_supg;
   
   % solving system
  
   U_0b(FreeDofs) = A_0b(FreeDofs,FreeDofs)\L_0b(FreeDofs);
   U_1b(FreeDofs) = A_1b(FreeDofs,FreeDofs)\L_1b(FreeDofs);
   U_1v(FreeDofs) = A_1v(FreeDofs,FreeDofs)\L_1v(FreeDofs);
   U_1m(FreeDofs) = A_1m(FreeDofs,FreeDofs)\L_1m(FreeDofs);
   
   U_2b(FreeDofs) = A_2b(FreeDofs,FreeDofs)\L_2b(FreeDofs);
   U_2v(FreeDofs) = A_2v(FreeDofs,FreeDofs)\L_2v(FreeDofs);
   U_2m(FreeDofs) = A_2m(FreeDofs,FreeDofs)\L_2m(FreeDofs);
   U_2e(FreeDofs) = A_2e(FreeDofs,FreeDofs)\L_2e(FreeDofs);
   
   U_supg(FreeDofs) = A_supg(FreeDofs,FreeDofs)\L_supg(FreeDofs);
   
   err(i,1) = L2Err_LFE(NewMesh,U_0b,P7O6(),d.U_EX_Handle,0,a);
   err(i,2) = L2Err_LFE(NewMesh,U_1b,P7O6(),d.U_EX_Handle,0,a);
   err(i,3) = L2Err_LFE(NewMesh,U_1v,P7O6(),d.U_EX_Handle,0,a);
   err(i,4) = L2Err_LFE(NewMesh,U_1m,P7O6(),d.U_EX_Handle,0,a);
   err(i,5) = L2Err_LFE(NewMesh,U_2b,P7O6(),d.U_EX_Handle,0,a);
   err(i,6) = L2Err_LFE(NewMesh,U_2v,P7O6(),d.U_EX_Handle,0,a);
   err(i,7) = L2Err_LFE(NewMesh,U_2m,P7O6(),d.U_EX_Handle,0,a);
   err(i,8) = L2Err_LFE(NewMesh,U_2e,P7O6(),d.U_EX_Handle,0,a);
   err(i,9) = L2Err_LFE(NewMesh,U_supg,P7O6(),d.U_EX_Handle,0,a);
   
   err2(i,1) = H1SErr_LFE(NewMesh,U_0b,P7O6(),d.GRAD_U_EX_Handle,0,a);
   err2(i,2) = H1SErr_LFE(NewMesh,U_1b,P7O6(),d.GRAD_U_EX_Handle,0,a);
   err2(i,3) = H1SErr_LFE(NewMesh,U_1v,P7O6(),d.GRAD_U_EX_Handle,0,a);
   err2(i,4) = H1SErr_LFE(NewMesh,U_1m,P7O6(),d.GRAD_U_EX_Handle,0,a);
   err2(i,5) = H1SErr_LFE(NewMesh,U_2b,P7O6(),d.GRAD_U_EX_Handle,0,a);
   err2(i,6) = H1SErr_LFE(NewMesh,U_2v,P7O6(),d.GRAD_U_EX_Handle,0,a);
   err2(i,7) = H1SErr_LFE(NewMesh,U_2m,P7O6(),d.GRAD_U_EX_Handle,0,a);
   err2(i,8) = H1SErr_LFE(NewMesh,U_2e,P7O6(),d.GRAD_U_EX_Handle,0,a);
   err2(i,9) = H1SErr_LFE(NewMesh,U_supg,P7O6(),d.GRAD_U_EX_Handle,0,a);
   
   Dofs(i)=size(NewMesh.Coordinates,1);
   
   %plot_LFE(U_up,NewMesh); colorbar;
   %plot_LFE(U_st,NewMesh);  colorbar;
   % plot_LFE(U_supg,NewMesh); colorbar;

   %plot_LFE(U_c,NewMesh); colorbar;
  end; 
  %U(FreeDofs)=M(FreeDofs,FreeDofs)\L_U(FreeDofs);
  %plot_LFE(U,NewMesh);  
  %plot_Mesh(NewMesh,'a')
  OFFSET = 0;
  tab=[1,2,3,4]; 
  
  xLim = [min(h)-OFFSET max(h)+OFFSET];
  
  yLim = [min(min(err))-OFFSET max(max(err))+OFFSET];
  
  yLim2 = [min(min(err2))-OFFSET max(max(err2))+OFFSET];
   
  fig = figure('Name','L^2 error 0 forms');
  plot(h,err(:,tab)); grid('on');
  set(gca,'XScale','log','YScale','log');
  xlabel('{\bf h}');
  ylabel('{\bf Error}');
  legend('(0,b)','(1,b)','(1,v)','(1,m)','Location','NorthWest');
  set(gca,'XLim',xLim,'YLim',yLim);
  Markers = '.ox+*sdv^<>ph';
  noMarkers = length(Markers);
  H1c = findobj(gca,'Type','line');
  linehandles = [H1c];
  for K = 1:length(linehandles)
     set(linehandles(K),'Marker',Markers(1+mod(K-1,noMarkers)));
  end
  add_Slope(gca,'South',1);
  add_Slope(gca,'Southeast',2);
  
  print -depsc './plots/L2rate0formsI1.eps'
  
  fig = figure('Name','H^1 semi error 0 forms');
  plot(h,err2(:,tab)); grid('on');
  set(gca,'XScale','log','YScale','log');
   Markers = '.ox+*sdv^<>ph';
  noMarkers = length(Markers);
  H1c = findobj(gca,'Type','line');
  linehandles = [H1c];
  for K = 1:length(linehandles)
     set(linehandles(K),'Marker',Markers(1+mod(K-1,noMarkers)));
  end
  xlabel('{\bf h}');
  ylabel('{\bf Error}');
  legend('(0,b)','(1,b)','(1,v)','(1,m)','Location','NorthWest'); 
  set(gca,'XLim',xLim,'YLim',yLim2);
  add_Slope(gca,'South',1);
  
  print -depsc './plots/H1semirate0formsI1.eps'

  tab=[5,6,7,8];
  fig = figure('Name','L^2 error 0 forms');
  plot(h,err(:,tab)); grid('on');
  set(gca,'XScale','log','YScale','log');
   Markers = '.ox+*sdv^<>ph';
  noMarkers = length(Markers);
  H1c = findobj(gca,'Type','line');
  linehandles = [H1c];
  for K = 1:length(linehandles)
     set(linehandles(K),'Marker',Markers(1+mod(K-1,noMarkers)));
  end
  xlabel('{\bf h}');
  ylabel('{\bf Error}');
  legend('(2,b)','(2,v)','(2,m)','(2,ex)','Location','NorthWest');
  set(gca,'XLim',xLim,'YLim',yLim);
  add_Slope(gca,'South',1);
  add_Slope(gca,'Southeast',2);
  
  print -depsc './plots/L2rate0formsI2.eps'
  
  fig = figure('Name','H^1 semi error 0 forms');
  plot(h,err2(:,tab)); grid('on');
  set(gca,'XScale','log','YScale','log');
   Markers = '.ox+*sdv^<>ph';
  noMarkers = length(Markers);
  H1c = findobj(gca,'Type','line');
  linehandles = [H1c];
  for K = 1:length(linehandles)
     set(linehandles(K),'Marker',Markers(1+mod(K-1,noMarkers)));
  end
  xlabel('{\bf h}');
  ylabel('{\bf Error}');
  legend('(2,b)','(2,v)','(2,m)','(2,ex)','Location','NorthWest')
  set(gca,'XLim',xLim,'YLim',yLim2);
  add_Slope(gca,'South',1);
  
  print -depsc './plots/H1semirate0formsI2.eps'
  
%   plot(h,l,'o-'); grid('on');
%   set(gca,'XScale','log','YScale','log');
%   xlabel('{\bf h}');
%   ylabel('{\bf Error}');
%   p = polyfit(log(h(NREFS-3:NREFS)),log(l(NREFS-3:NREFS,1)),1);
%   add_Slope(gca,'Southwest',p(1),'g-');
%   p = polyfit(log(h(NREFS-3:NREFS)),log(l(NREFS-3:NREFS,2)),1);
%   add_Slope(gca,'Southeast',p(1),'g-');
  
clear all;
