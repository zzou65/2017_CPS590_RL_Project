 % Initialize constants
  clear Mesh
 
  JIG=1;
  d1=0.1;                    %impact of supg modification
  d2=0.1;
  a=10^0                %amount of diffusivity
  NREFS =6;
  % select test code
  d=getData(8);
    
  QuadRule = P7O6(); 
  
  Mesh.Coordinates =d.Coordinates; 
  Mesh.Elements = d.Elements;
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = d.boundtype;
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  
  err=zeros(NREFS,6);
  h=zeros(NREFS,1);
  Dofs=zeros(NREFS,1);
  
  for i = 1:NREFS
      i
   
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
   [qr_bnd qr_inn]=split_QuadRule(P10O4());
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

   weights3=[0.05 0.05 0.05 4/30 4/30 4/30, 0.45];
   [M3,UP3,C3]=assemMat_UpQFE(New_Mesh,d.V_Handle,weights3);
   
   % standard
   B=assemMat_QFE(New_Mesh,@STIMA_Conv_QFE,d.V_Handle,P3O3());
   
   A_u2=a*A+(M3*UP3+C3);
   A_u3=a*A+B3;
   A_u4=a*A+B4;
   A_u5=a*A+B5; 
   A_s=a*A+B;
   
   L = assemLoad_QFE(New_Mesh,P3O3(),d.SOL_Handle,a);
   
   %Incorporate Dirichlet boundary data
   
   [U_s,FreeDofs] = assemDir_QFE(New_Mesh,[-1 -2],d.U_EX_Handle,a);
   U_u2=U_s;
   U_u3=U_s;
   U_u4=U_s;
   U_u5=U_s;
    
   L_u2 = L - A_u2*U_u2;
   L_u3 = L - A_u3*U_u3;
   L_u4 = L - A_u4*U_u4;
   L_u5 = L - A_u5*U_u5;
   L_s = L - A_s*U_s;
   
   % Solve the linear system
 
   U_u2(FreeDofs) = A_u2(FreeDofs,FreeDofs)\L_u2(FreeDofs);
   U_u3(FreeDofs) = A_u3(FreeDofs,FreeDofs)\L_u3(FreeDofs);
   U_u4(FreeDofs) = A_u4(FreeDofs,FreeDofs)\L_u4(FreeDofs);
   U_u5(FreeDofs) = A_u5(FreeDofs,FreeDofs)\L_u5(FreeDofs);
   U_s(FreeDofs) = A_s(FreeDofs,FreeDofs)\L_s(FreeDofs);
   q=Duffy(TProd(gauleg(0,1,10)));
   err(i,2) = L2Err_QFE(New_Mesh,U_u2,q,d.U_EX_Handle,0,a);
   err(i,3) = L2Err_QFE(New_Mesh,U_u3,q,d.U_EX_Handle,0,a);
   err(i,4) = L2Err_QFE(New_Mesh,U_u4,P7O6(),d.U_EX_Handle,0,a);
   err(i,5) = L2Err_QFE(New_Mesh,U_u5,P7O6(),d.U_EX_Handle,0,a);
   err(i,6) = L2Err_QFE(New_Mesh,U_s,P7O6(),d.U_EX_Handle,0,a);
   
   Dofs(i)=size(New_Mesh.Coordinates,1);  
  
%     plot_QFE(U_u4,New_Mesh); colorbar
%     plot_QFE(U_u1,New_Mesh); colorbar
%        plot_QFE(U_df,New_Mesh); colorbar
%     plot_QFE(U_u3,New_Mesh); colorbar
%       plotLine_QFE(U_df,New_Mesh,[0,1], [1,0]);
%      plotLine_QFE(U_u3,New_Mesh,[0,1], [1,0]);
%      plotLine_QFE(U_s,New_Mesh,[0,1], [1,0]);
%      plot_QFE(U_s,New_Mesh); colorbar
  %   plot_QFE(U_supg,New_Mesh); colorbar
 end; 
  h=Dofs;
  fig = figure('Name','Discretization error');
  plot(h,err(:,2),'ro-',h,err(:,3),'go-',h,err(:,4),'yo-',h,err(:,5),'co-',h,err(:,6),'mo-'); grid('on');
  set(gca,'XScale','log','YScale','log');
  xlabel('{\bf h}');
  ylabel('{\bf Error}');
  legend('P704','P704','NCC(5)','P6O4','FEM','Location','NorthEast');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,3)),1);
  add_Slope(gca,'East',p(1),'r-');  
  p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,3)),1);
  add_Slope(gca,'Southwest',p(1),'g-');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,4)),1);
  add_Slope(gca,'West',p(1),'y-');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,5)),1);
  add_Slope(gca,'NorthWest',p(1),'c-');
  p = polyfit(log(h(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,6)),1);
  add_Slope(gca,'NorthEast',p(1),'m-');
clear all;
