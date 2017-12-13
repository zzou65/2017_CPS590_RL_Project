% Initialize constants
  clear Mesh
 
  JIG=1;
  d1=0.1;                    %impact of supg modification
  d2=0.1;
  a=10^0                %amount of diffusivity
  NREFS =4;
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
  
  err=zeros(NREFS,5);
  h=zeros(NREFS,1);
  Dofs=zeros(NREFS,1);
  
  for i = 1:NREFS
   
   %refine Mesh   
   Mesh = refine_REG(Mesh);
   Mesh=add_Edge2Elem(Mesh);
        
   % Mesh preprocessing
   switch(JIG)
        case 1
          New_Mesh = Mesh;      
        case 2
          Loc = get_BdEdges(Mesh);
          Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
          FixedPos = zeros(size(Mesh.Coordinates,1),1);
          FixedPos(Loc) = 1;
          New_Mesh = jiggle(Mesh,FixedPos);   
        case 3
          Loc = get_BdEdges(Mesh);
          Loc = unique([Mesh.Edges(Loc,1); Mesh.Edges(Loc,2)]);
          FixedPos = zeros(size(Mesh.Coordinates,1),1);
          FixedPos(Loc) = 1;
          New_Mesh = smooth(Mesh,FixedPos);
   end
        
   h(i)=get_MeshWidth(Mesh);
   
   %Laplace
   %A = assemMat_QFE(New_Mesh,@STIMA_Lapl_QFE_quad,P3O3());

   A = assemMat_QFE(New_Mesh,@STIMA_Lapl_QFE);

   %    weights1=[0 0 0 1/3 1/3 1/3 0];
   %    weights2=[1/12 1/12 1/12 0 0 0 3/4];
   %  %  weights2=[1/12 1/12 1/12 1/4 1/4 1/4 0];
   %    weights3=[0.05 0.05 0.05 4/30 4/30 4/30, 0.45];
   %    weights4=[1/3 1/3 1/3 0 0 0 0];
   %    [M1,UP1,C1]=assemMat_UpQFE(New_Mesh,d.V_Handle,weights1);
   %    [M2,UP2,C2]=assemMat_UpQFE(New_Mesh,d.V_Handle,weights2);
   %    [M3,UP3,C3]=assemMat_UpQFE(New_Mesh,d.V_Handle,weights3);
   %    [M4,UP4,C4]=assemMat_UpQFE(New_Mesh,d.V_Handle,weights4);
   %    A_u1=a*A+(M1*UP1);
   %    A_u2=a*A+(M2*UP2+C2);
   %    A_u3=a*A+(M3*UP3+C3);
   %   A_u4=a*A+(M4*UP4+C4);
   % new version
   % weight1
   
   [qr_bnd qr_inn]=split_QuadRule(P3O3());
   if (~isempty(qr_inn.x))
        B1=assemMat_QFE(New_Mesh,@STIMA_Conv_QFE,d.V_Handle,qr_inn);
        if (~isempty(qr_bnd.x))
            B1=B1+assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   else
        if (~isempty(qr_bnd.x))
            B1=assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   end
   
   
   % weight2
   [qr_bnd qr_inn]=split_QuadRule(P4O3());
    if (~isempty(qr_inn.x))
        B2=assemMat_QFE(New_Mesh,@STIMA_Conv_QFE,d.V_Handle,qr_inn);
        if (~isempty(qr_bnd.x))
            B2=B2+assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   else
        if (~isempty(qr_bnd.x))
            B2=assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   end
   
   % weight3
   [qr_bnd qr_inn]=split_QuadRule(P7O4());
    if (~isempty(qr_inn.x))
        B3=assemMat_QFE(New_Mesh,@STIMA_Conv_QFE,d.V_Handle,qr_inn);
        if (~isempty(qr_bnd.x))
            B3=B3+assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   else
        if (~isempty(qr_bnd.x))
            B3=assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   end
   
   % weight4
%   [qr_bnd qr_inn]=split_QuadRule(Duffy(TProd(gaulob(0,1,5))));
   [qr_bnd qr_inn]=split_QuadRule(NCC(5));
   if (~isempty(qr_inn.x))
        B4=assemMat_QFE(New_Mesh,@STIMA_Conv_QFE,d.V_Handle,qr_inn);
        if (~isempty(qr_bnd.x))
            B4=B4+assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   else
        if (~isempty(qr_bnd.x))
            B4=assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   end
   
   % weight5
   [qr_bnd qr_inn]=split_QuadRule(P6O4());
   if (~isempty(qr_inn.x))
        B5=assemMat_QFE(New_Mesh,@STIMA_Conv_QFE,d.V_Handle,qr_inn);
        if (~isempty(qr_bnd.x))
            B5=B5+assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   else
        if (~isempty(qr_bnd.x))
            B5=assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   end
   
   % weight6
   [qr_bnd qr_inn]=split_QuadRule(P10O5());
   if (~isempty(qr_inn.x))
        B6=assemMat_QFE(New_Mesh,@STIMA_Conv_QFE,d.V_Handle,qr_inn);
        if (~isempty(qr_bnd.x))
            B6=B6+assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   else
        if (~isempty(qr_bnd.x))
            B6=assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   end
   A_u1=a*A+B1;
   A_u2=a*A+B2;
   A_u3=a*A+B3;
   A_u4=a*A+B4;
   A_u5=a*A+B5;
   A_u6=a*A+B6;

   B_supg=assemMat_QFE(New_Mesh,@STIMA_SUPG_QFE,P3O3(), d.V_Handle,a,d1,d2);
   B=assemMat_QFE(New_Mesh,@STIMA_Conv_QFE,d.V_Handle,P3O3());
  
   B_df=assemMat_QFE(New_Mesh,@STIMA_UPQFE_DF,d.V_Handle);
   M_df=assemMat_QFE(New_Mesh,@MASS_QFEquad,P3O3());
   B_df=M_df*B_df;
   
%    B_df2=assemMat_QFELFE(New_Mesh,@STIMA_UPQFE_DF,d.V_Handle);
%    M_df=assemMat_QFE(New_Mesh,@MASS_QFEquad,P7O4());
%    B_df=M_df*B_df;
   
   A_supg=a*A+B+B_supg;
   A_s=a*A+B;
   A_df=a*A+B_df;
   L = assemLoad_QFE(New_Mesh,P3O3(),d.SOL_Handle,a);
   L_supg=assemLoad_QFE_SUPG(New_Mesh,P3O3(),d.V_Handle,d.SOL_Handle,a,d1,d2);
   
   %Incorporate Dirichlet boundary data
   
   [U_s,FreeDofs] = assemDir_QFE(New_Mesh,[-1 -2],d.U_EX_Handle,a);
   U_u1=U_s;
   U_u2=U_s;
   U_u3=U_s;
   U_u4=U_s;
   U_u5=U_s;
   U_u6=U_s;
   U_supg=U_s;
   U_df=U_s;
    
   L_u1 = L - A_u1*U_u1;
   L_u2 = L - A_u2*U_u2;
   L_u3 = L - A_u3*U_u3;
   L_u4 = L - A_u4*U_u4;
   L_u5 = L - A_u5*U_u5;
   L_u6 = L - A_u6*U_u6;
   L_s = L - A_s*U_s;
   L_df = L - A_df*U_df;
   L_supg = L+L_supg - A_supg*U_supg;
   
   % Solve the linear system
 
   U_u1(FreeDofs) = A_u1(FreeDofs,FreeDofs)\L_u1(FreeDofs);
   U_u2(FreeDofs) = A_u2(FreeDofs,FreeDofs)\L_u2(FreeDofs);
   U_u3(FreeDofs) = A_u3(FreeDofs,FreeDofs)\L_u3(FreeDofs);
   U_u4(FreeDofs) = A_u4(FreeDofs,FreeDofs)\L_u4(FreeDofs);
   U_u5(FreeDofs) = A_u5(FreeDofs,FreeDofs)\L_u5(FreeDofs);
   U_u6(FreeDofs) = A_u6(FreeDofs,FreeDofs)\L_u6(FreeDofs);
   U_s(FreeDofs) = A_s(FreeDofs,FreeDofs)\L_s(FreeDofs);
   U_df(FreeDofs) = A_df(FreeDofs,FreeDofs)\L_df(FreeDofs);
   U_supg(FreeDofs) = A_supg(FreeDofs,FreeDofs)\L_supg(FreeDofs);
   %   
   err(i,1) = L2Err_QFE(New_Mesh,U_u1,P7O6(),d.U_EX_Handle,0,a);
   err(i,2) = L2Err_QFE(New_Mesh,U_u2,P7O6(),d.U_EX_Handle,0,a);
   err(i,3) = L2Err_QFE(New_Mesh,U_u3,P7O6(),d.U_EX_Handle,0,a);
   err(i,4) = L2Err_QFE(New_Mesh,U_u4,P7O6(),d.U_EX_Handle,0,a);
   err(i,5) = L2Err_QFE(New_Mesh,U_u5,P7O6(),d.U_EX_Handle,0,a);
   err(i,6) = L2Err_QFE(New_Mesh,U_df,P7O6(),d.U_EX_Handle,0,a);
%   err(i,5) = L2Err_QFE(New_Mesh,U_s,P7O6(),d.U_EX_Handle,0,a);
%   err(i,6) = L2Err_QFE(New_Mesh,U_supg,P7O6(),d.U_EX_Handle,0,a)
   err2(i,1) = H1SErr_QFE(New_Mesh,U_u1,P7O6(),d.GRAD_U_EX_Handle,0,a);
   err2(i,2) = H1SErr_QFE(New_Mesh,U_u2,P7O6(),d.GRAD_U_EX_Handle,0,a);
   err2(i,3) = H1SErr_QFE(New_Mesh,U_u3,P7O6(),d.GRAD_U_EX_Handle,0,a);
   err2(i,4) = H1SErr_QFE(New_Mesh,U_u4,P7O6(),d.GRAD_U_EX_Handle,0,a);
   err2(i,5) = H1SErr_QFE(New_Mesh,U_u5,P7O6(),d.GRAD_U_EX_Handle,0,a);
   err2(i,6) = H1SErr_QFE(New_Mesh,U_df,P7O6(),d.GRAD_U_EX_Handle,0,a);
   Dofs(i)=size(New_Mesh.Coordinates,1);  
  
%     plot_QFE(U_u4,New_Mesh); colorbar
%     plot_QFE(U_u1,New_Mesh); colorbar
%        plot_QFE(U_df,New_Mesh); colorbar
%     plot_QFE(U_u3,New_Mesh); colorbar
%       plotLine_QFE(U_df,New_Mesh,[0,1], [1,0]);
%      plotLine_QFE(U_u3,New_Mesh,[0,1], [1,0]);
   %  plotLine_QFE(U_s,New_Mesh,[0,1], [1,0]);
   %  plot_QFE(U_s,New_Mesh); colorbar
  %   plot_QFE(U_supg,New_Mesh); colorbar
 end; 
  fig = figure('Name','Discretization error');
  plot(h,err(:,4),'go-',h,err(:,1),'ro-',h,err(:,2),'bo-',h,err(:,3),'yo-',h,err(:,5),'co-',h,err(:,6),'mo-'); grid('on');
  set(gca,'XScale','log','YScale','log');
  xlabel('{\bf h}');
  ylabel('{\bf Error}');
  legend('L^2-error u (up1)','L^2-error u (up2)','L^2-error u (up3)',...
      'L^2-error u (up4)','L^2-error u','L^2-error u(supg)','Location','NorthEast');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,4)),1);
  add_Slope(gca,'West',p(1),'g-');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,1)),1);
  add_Slope(gca,'SouthEast',p(1),'r-');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,2)),1);
  add_Slope(gca,'East',p(1),'b-');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,3)),1);
  add_Slope(gca,'Southwest',p(1),'y-');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,5)),1);
  add_Slope(gca,'NorthWest',p(1),'c-');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,6)),1);
  add_Slope(gca,'NorthEast',p(1),'m-');
  
  fig = figure('Name','Discretization error');
  plot(h,err2(:,4),'go-',h,err2(:,1),'ro-',h,err2(:,2),'bo-',h,err2(:,3),'yo-',h,err2(:,5),'co-',h,err2(:,6),'mo-'); grid('on');
  set(gca,'XScale','log','YScale','log');
  xlabel('{\bf h}');
  ylabel('{\bf Error}');
  legend('H^1-error u (up1)','H^1-error u (up2)','H^1-error u (up3)',...
      'H^1-error u (up4)','H^1-error u','H^1-error u(supg)','Location','NorthEast');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err2(NREFS-3:NREFS,4)),1);
  add_Slope(gca,'West',p(1),'g-');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err2(NREFS-3:NREFS,1)),1);
  add_Slope(gca,'SouthEast',p(1),'r-');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err2(NREFS-3:NREFS,2)),1);
  add_Slope(gca,'East',p(1),'b-');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err2(NREFS-3:NREFS,3)),1);
  add_Slope(gca,'Southwest',p(1),'y-');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err2(NREFS-3:NREFS,5)),1);
  add_Slope(gca,'NorthWest',p(1),'c-');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err2(NREFS-3:NREFS,6)),1);
  add_Slope(gca,'NorthEast',p(1),'m-');
clear all;
