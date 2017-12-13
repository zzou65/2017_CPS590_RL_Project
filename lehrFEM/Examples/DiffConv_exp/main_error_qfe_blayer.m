% driver routine for several diffusion convection discretizations
% see below
   clear Mesh
  % Initialize constants
  JIG=2;
  d1=1;                      %impact of supg modification
  d2=1;
  a=10^0;               %amount of diffusivity
  NREFS =5;
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
   A = assemMat_QFE(NewMesh,@STIMA_Lapl_QFE);
   
   % linear interpolation
%    M_QFELFE_b=assemMat_MASS_QFELFE(NewMesh,P1O2());
%    M_QFELFE_v=assemMat_MASS_QFELFE(NewMesh,P3O2());
   M_QFELFE_m=assemMat_MASS_QFELFE(NewMesh,P3O3());
   M_QFELFE_e=assemMat_MASS_QFELFE(NewMesh,P4O3());
   
   I_1=assemMat_QFELFE(NewMesh, @STIMA_UPQFE_LFE, d.V_Handle);
   
%    D_1b=M_QFELFE_b'*I_1;
%    D_1v=M_QFELFE_v'*I_1;
   D_1m=M_QFELFE_m'*I_1;
   D_1e=M_QFELFE_e'*I_1;

   % standard Galerkin D_0b=D0c
   D_0b=assemMat_QFE(NewMesh,@STIMA_Conv_QFE ,d.V_Handle,P3O3());
    
   % SUPG
%   D_supg2 = assemMat_QFE(NewMesh,@STIMA_SUPG_QFE,P3O2(), d.V_Handle,a,d1,d2);
%   D_supg  = D_0b+D_supg2;
   
   % quadratic interpolation
%   M_QFE_b=assemMat_QFE(NewMesh,@MASS_QFEquad,P1O2());
%   M_QFE_v=assemMat_QFE(NewMesh,@MASS_QFEquad,P3O2());
   M_QFE_vmb=assemMat_QFE(NewMesh,@MASS_QFEquad,P7O4());
   M_QFE_e=assemMat_QFE(NewMesh,@MASS_QFEquad,P7O6());
   
   I_2=assemMat_QFE(NewMesh,@STIMA_UPQFE_QFE,d.V_Handle);
   %I_2=assemMat_QFE(NewMesh,@STIMA_UPQFE_QFE,d.V_Handle,h);
   
%    D_2b=M_QFE_b*I_2;
%    D_2v=M_QFE_v*I_2;
   D_2vmb=M_QFE_vmb*I_2;
   D_2e=M_QFE_e*I_2;
   
   
   [qr_bnd qr_inn]=split_QuadRule(P1O2());
   if (~isempty(qr_inn.x))
        B1=assemMat_QFE(NewMesh,@STIMA_Conv_QFE,d.V_Handle,qr_inn);
        if (~isempty(qr_bnd.x))
            B1=B1+assemMat_UPQFE2(NewMesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   else
        if (~isempty(qr_bnd.x))
            B1=assemMat_UPQFE2(NewMesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   end
   
   
   % weight2
   [qr_bnd qr_inn]=split_QuadRule(P6O4());
    if (~isempty(qr_inn.x))
        B2=assemMat_QFE(NewMesh,@STIMA_Conv_QFE,d.V_Handle,qr_inn);
        if (~isempty(qr_bnd.x))
            B2=B2+assemMat_UPQFE2(NewMesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   else
        if (~isempty(qr_bnd.x))
            B2=assemMat_UPQFE2(NewMesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   end
   
   % weight3
   [qr_bnd qr_inn]=split_QuadRule(P7O4());
    if (~isempty(qr_inn.x))
        B3=assemMat_QFE(NewMesh,@STIMA_Conv_QFE,d.V_Handle,qr_inn);
        if (~isempty(qr_bnd.x))
            B3=B3+assemMat_UPQFE2(NewMesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   else
        if (~isempty(qr_bnd.x))
            B3=assemMat_UPQFE2(NewMesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   end

   %Laplace + convection
   A1=a*A+B1;
   A2=a*A+B2;
   A3=a*A+B3;
   A_0b =(a*A+D_0b);
%    A_1b =(a*A+D_1b);
%    A_1v =(a*A+D_1v);
   A_1m =(a*A+D_1m);
   A_1e =(a*A+D_1e);   
%    A_2b =(a*A+D_2b);
%    A_2v =(a*A+D_2v);
   A_2vmb =(a*A+D_2vmb);
   A_2e =(a*A+D_2e);
%   A_supg =(a*A+D_supg);
   % source term
  
   L = assemLoad_QFE(NewMesh,P7O6(),d.SOL_Handle,a);
   L_U = assemLoad_QFE(NewMesh,P7O6(),d.U_EX_Handle,a);
   L_supg=assemLoad_QFE_SUPG(NewMesh,P7O6(),d.V_Handle,d.SOL_Handle,a,d1,d2);
   
   % Direchlet boundary
  
   [U_0b,FreeDofs] = assemDir_QFE(NewMesh,-1,d.U_EX_Handle,a);
%    U_1b=U_0b;
%    U_1v=U_0b;
   U_1m=U_0b;
   U_1e=U_0b;
%    U_2b=U_0b;
%    U_2v=U_0b;
   U_2vmb=U_0b;
   U_2e=U_0b;
%    U_supg=U_0b;
   U1=U_0b;
   U2=U_0b;
   U3=U_0b;
   
   L_0b = L - A_0b*U_0b;
%    L_1b = L - A_1b*U_1b;
%    L_1v = L - A_1v*U_1v;
   L_1m = L - A_1m*U_1m;
   L_1e = L - A_1e*U_1e;
%    L_2b = L - A_2b*U_2b;
%    L_2v = L - A_2v*U_2v;
   L_2vmb = L - A_2vmb*U_2vmb;
   L_2e = L - A_2e*U_2e;
   
   L1 = L - A1*U1;
   L2 = L - A2*U2;
   L3 = L - A3*U3;

%    L_supg = (L+L_supg) - A_supg*U_supg;
   
   % solving system
  
   U_0b(FreeDofs) = A_0b(FreeDofs,FreeDofs)\L_0b(FreeDofs);
%    U_1b(FreeDofs) = A_1b(FreeDofs,FreeDofs)\L_1b(FreeDofs);
%    U_1v(FreeDofs) = A_1v(FreeDofs,FreeDofs)\L_1v(FreeDofs);
   U_1m(FreeDofs) = A_1m(FreeDofs,FreeDofs)\L_1m(FreeDofs);
   U_1e(FreeDofs) = A_1e(FreeDofs,FreeDofs)\L_1e(FreeDofs);
   
%    U_2b(FreeDofs) = A_2b(FreeDofs,FreeDofs)\L_2b(FreeDofs);
%    U_2v(FreeDofs) = A_2v(FreeDofs,FreeDofs)\L_2v(FreeDofs);
   U_2vmb(FreeDofs) = A_2vmb(FreeDofs,FreeDofs)\L_2vmb(FreeDofs);
   U_2e(FreeDofs) = A_2e(FreeDofs,FreeDofs)\L_2e(FreeDofs);
   
   U1(FreeDofs) = A1(FreeDofs,FreeDofs)\L1(FreeDofs);
   U2(FreeDofs) = A2(FreeDofs,FreeDofs)\L2(FreeDofs);
   U3(FreeDofs) = A3(FreeDofs,FreeDofs)\L3(FreeDofs);
   
%    U_supg(FreeDofs) = A_supg(FreeDofs,FreeDofs)\L_supg(FreeDofs);
   
   err(i,1) = L2Err_QFE(NewMesh,U_0b,P7O6(),d.U_EX_Handle,0,a);
   %err(i,2) = L2Err_QFE(NewMesh,U_1b,P7O6(),d.U_EX_Handle,0,a);
   %err(i,3) = L2Err_QFE(NewMesh,U_1v,P7O6(),d.U_EX_Handle,0,a);
   err(i,2) = L2Err_QFE(NewMesh,U_1m,P7O6(),d.U_EX_Handle,0,a);
   err(i,3) = L2Err_QFE(NewMesh,U_1e,P7O6(),d.U_EX_Handle,0,a);
%    err(i,5) = L2Err_QFE(NewMesh,U_2b,P7O6(),d.U_EX_Handle,0,a);
%    err(i,6) = L2Err_QFE(NewMesh,U_2v,P7O6(),d.U_EX_Handle,0,a);
   err(i,4) = L2Err_QFE(NewMesh,U_2vmb,P7O6(),d.U_EX_Handle,0,a);
   err(i,5) = L2Err_QFE(NewMesh,U_2e,P7O6(),d.U_EX_Handle,0,a);
%   err(i,9) = L2Err_QFE(NewMesh,U_supg,P7O6(),d.U_EX_Handle,0,a);
   err(i,6) = L2Err_QFE(NewMesh,U1,P7O6(),d.U_EX_Handle,0,a);
   err(i,7) = L2Err_QFE(NewMesh,U2,P7O6(),d.U_EX_Handle,0,a);
   err(i,8) = L2Err_QFE(NewMesh,U3,P7O6(),d.U_EX_Handle,0,a);
   
   err2(i,1) = L2Err_QFEmod(NewMesh,U_0b,P7O6(),d.U_EX_Handle,0,a);
  % err2(i,2) = H1SErr_QFE(NewMesh,U_1b,P7O6(),d.GRAD_U_EX_Handle,0,a);
  % err2(i,3) = H1SErr_QFE(NewMesh,U_1v,P7O6(),d.GRAD_U_EX_Handle,0,a);
   err2(i,2) = L2Err_QFEmod(NewMesh,U_1m,P7O6(),d.U_EX_Handle,0,a);
   err2(i,3) = L2Err_QFEmod(NewMesh,U_1e,P7O6(),d.U_EX_Handle,0,a);
   %err2(i,5) = H1SErr_QFE(NewMesh,U_2b,P7O6(),d.GRAD_U_EX_Handle,0,a);
   %err2(i,6) = H1SErr_QFE(NewMesh,U_2v,P7O6(),d.GRAD_U_EX_Handle,0,a);
   err2(i,4) = H1SErr_QFE(NewMesh,U_2vmb,P7O6(),d.GRAD_U_EX_Handle,0,a);
   err2(i,5) = L2Err_QFEmod(NewMesh,U_2e,P7O6(),d.U_EX_Handle,0,a);
   %err2(i,9) = H1SErr_QFE(NewMesh,U_supg,P7O6(),d.GRAD_U_EX_Handle,0,a);
   err2(i,6) = L2Err_QFEmod(NewMesh,U1,P7O6(),d.U_EX_Handle,0,a);
   err2(i,7) = L2Err_QFEmod(NewMesh,U2,P7O6(),d.U_EX_Handle,0,a);
   err2(i,8) = L2Err_QFEmod(NewMesh,U3,P7O6(),d.U_EX_Handle,0,a);
   
   Dofs(i)=size(NewMesh.Coordinates,1);
   
   %plot_LFE(U_up,NewMesh); colorbar;
   plot_QFE(U_2e, NewMesh);  colorbar;
   plot_QFE(U2, NewMesh); colorbar;
   plot_QFE(U3, NewMesh); colorbar;

   %plot_LFE(U_c,NewMesh); colorbar;
  end; 
  %U(FreeDofs)=M(FreeDofs,FreeDofs)\L_U(FreeDofs);
  %plot_LFE(U,NewMesh);  
  %plot_Mesh(NewMesh,'a')
  OFFSET = 0;
  tab=1:8; 
  
  xLim = [min(h)-OFFSET max(h)+OFFSET];
  
  yLim = [min(min(err))-OFFSET max(max(err))+OFFSET];
  
  yLim2 = [min(min(err2))-OFFSET max(max(err2))+OFFSET];
   
  fig = figure('Name','L^2 error 0 forms');
  plot(h,err(:,tab),'o-'); grid('on');
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
 legend('stand','(1,m)','(1,vb)','(2,vmb)','(2,ex)','(2,m)','(1b,vb)','(2b,vmb)','Location','NorthWest');
  set(gca,'XLim',xLim,'YLim',yLim);
  add_Slope(gca,'South',2);
  add_Slope(gca,'Southeast',3);
%   p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,4)),1);
%   add_Slope(gca,'Southwest',p(1),'g-');
  
  print -depsc './plots/L2rate0forms_qfe_blayer.eps'
  
  fig = figure('Name','error 0 forms excl blayer');
  plot(h,err2(:,tab),'o-'); grid('on');
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
legend('stand','(1,m)','(1,vb)','(2,vmb)','(2,ex)','(2,m)','(1b,vb)','(2b,vmb)','Location','NorthWest');  %set(gca,'XLim',xLim,'YLim',yLim2);
  add_Slope(gca,'east',2);
%   p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,4)),1);
%   add_Slope(gca,'Southwest',p(1),'g-');
  
  print -depsc './plots/L2modrate0forms_qfe_blayer.eps'

%   tab=[5,6,7,8,9];
%   fig = figure('Name','L^2 error 0 forms');
%   plot(h,err(:,tab)); grid('on');
%   set(gca,'XScale','log','YScale','log');
%   xlabel('{\bf h}');
%   ylabel('{\bf Error}');
%   legend('(2,b)','(2,v)','(2,m)','(2,e)','(SUPG)','Location','NorthWest');
%   set(gca,'XLim',xLim,'YLim',yLim);
%   add_Slope(gca,'South',2);
%   add_Slope(gca,'Southeast',3);
% %   p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,7)),1);
% %   add_Slope(gca,'Southwest',p(1),'g-');
%   
%   print -deps 'L2rate0formsI2_qfe.eps'
%   
%   fig = figure('Name','H^1 semi error 0 forms');
%   plot(h,err2(:,tab)); grid('on');
%   set(gca,'XScale','log','YScale','log');
%   xlabel('{\bf h}');
%   ylabel('{\bf Error}');
%   legend('(2,b)','(2,v)','(2,m)','(2,e)','(SUPG)','Location','NorthWest')
%   set(gca,'XLim',xLim,'YLim',yLim2);
%   add_Slope(gca,'South',2);
% %   p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,7)),1);
% %   add_Slope(gca,'Southwest',p(1),'g-');
%   
%   print -deps 'H1semirate0formsI2_qfe.eps'
%   
clear all;