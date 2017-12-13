 % Initialize constants
  clear Mesh
 
  JIG=1;
  d1=1;                    %impact of supg modification
  d2=1;
  a=10^-4               %amount of diffusivity
  NREFS =4;
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
  
  err=zeros(NREFS,5);
  h=zeros(NREFS,1);
  Dofs=zeros(NREFS,1);
  
  Mesh = refine_REG(Mesh);
  Mesh=add_Edge2Elem(Mesh);
  
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
   
   % Laplace
   
   A = assemMat_QFE(New_Mesh,@STIMA_Lapl_QFE);
   
   % Mass Matrix
   
   M = assemMat_QFE(New_Mesh,@MASS_QFE);
   
   % Convection
   
%    weights1=[0 0 0 1/3 1/3 1/3 0];
%    weights2=[1/12 1/12 1/12 0 0 0 3/4];
%    %  weights2=[1/12 1/12 1/12 1/4 1/4 1/4 0];
%    weights3=[0.05 0.05 0.05 4/30 4/30 4/30, 0.45];
%    weights4=[1/3 1/3 1/3 0 0 0 0];
%    
%    [M1,UP1,C1]=assemMat_UpQFE(New_Mesh,V_Handle,weights1);
%    [M2,UP2,C2]=assemMat_UpQFE(New_Mesh,V_Handle,weights2);
%    [M3,UP3,C3]=assemMat_UpQFE(New_Mesh,V_Handle,weights3);
%    [M4,UP4,C4]=assemMat_UpQFE(New_Mesh,V_Handle,weights4);
%    A_u1=a*A+(M1*UP1);
%    A_u2=a*A+(M2*UP2+C2);
%    A_u3=a*A+(M3*UP3+C3);
%    A_u4=a*A+(M4*UP4+C4);
   [qr_bnd qr_inn]=split_QuadRule(P3O3());
   if (~isempty(qr_inn.x))
        B1=assemMat_QFE(New_Mesh,@STIMA_Conv_QFE,V_Handle,qr_inn);
        if (~isempty(qr_bnd.x))
            B1=B1+assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,V_Handle,qr_bnd);
        end
   else
        if (~isempty(qr_bnd.x))
            B1=assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,V_Handle,qr_bnd);
        end
   end
   
   
   % weight2
   [qr_bnd qr_inn]=split_QuadRule(P4O3());
    if (~isempty(qr_inn.x))
        B2=assemMat_QFE(New_Mesh,@STIMA_Conv_QFE,V_Handle,qr_inn);
        if (~isempty(qr_bnd.x))
            B2=B2+assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,V_Handle,qr_bnd);
        end
   else
        if (~isempty(qr_bnd.x))
            B2=assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,V_Handle,qr_bnd);
        end
   end
   
   % weight3
   [qr_bnd qr_inn]=split_QuadRule(P7O4());
    if (~isempty(qr_inn.x))
        B3=assemMat_QFE(New_Mesh,@STIMA_Conv_QFE,V_Handle,qr_inn);
        if (~isempty(qr_bnd.x))
            B3=B3+assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,V_Handle,qr_bnd);
        end
   else
        if (~isempty(qr_bnd.x))
            B3=assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,V_Handle,qr_bnd);
        end
   end
   
   % weight4
   %[qr_bnd qr_inn]=split_QuadRule(Duffy(TProd(gaulob(0,1,3))));
   [qr_bnd qr_inn]=split_QuadRule(P3O2());
   if (~isempty(qr_inn.x))
        B4=assemMat_QFE(New_Mesh,@STIMA_Conv_QFE,V_Handle,qr_inn);
        if (~isempty(qr_bnd.x))
            B4=B4+assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,V_Handle,qr_bnd);
        end
   else
        if (~isempty(qr_bnd.x))
            B4=assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,V_Handle,qr_bnd);
        end
   end
   
   % weight5
   [qr_bnd qr_inn]=split_QuadRule(P6O4());
   if (~isempty(qr_inn.x))
        B5=assemMat_QFE(New_Mesh,@STIMA_Conv_QFE,V_Handle,qr_inn);
        if (~isempty(qr_bnd.x))
            B5=B5+assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,V_Handle,qr_bnd);
        end
   else
        if (~isempty(qr_bnd.x))
            B5=assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,V_Handle,qr_bnd);
        end
   end
   
   % weight6
   [qr_bnd qr_inn]=split_QuadRule(P10O5());
   if (~isempty(qr_inn.x))
        B6=assemMat_QFE(New_Mesh,@STIMA_Conv_QFE,V_Handle,qr_inn);
        if (~isempty(qr_bnd.x))
            B6=B6+assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,V_Handle,qr_bnd);
        end
   else
        if (~isempty(qr_bnd.x))
            B6=assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,V_Handle,qr_bnd);
        end
   end
   A_u1=a*A+B1;
   A_u2=a*A+B2;
   A_u3=a*A+B3;
   A_u4=a*A+B4;
   A_u5=a*A+B5;
   A_u6=a*A+B6;

   B_supg=assemMat_QFE(New_Mesh,@STIMA_SUPG_QFE,P3O3(), V_Handle,a,d1,d2);
   B=assemMat_QFE(New_Mesh,@STIMA_Conv_QFE,V_Handle,P3O3());
   
   B_df=assemMat_QFE(New_Mesh,@STIMA_UPQFE_DF,V_Handle);
   M_df=assemMat_QFE(New_Mesh,@MASS_QFEquad,P3O2());
   B_df=M_df*B_df;
   
   A_supg=a*A+B+B_supg;
   A_s=a*A+B;
   A_df=a*A+B_df;

   % discrete inf-sup-condition
   
   FreeDofs = FreeDofs_QFE(New_Mesh,[-1,-2]);
   
   H=assemMat_QFE_BdLayer(New_Mesh,@STIMA_InfSup_QFE,QuadRule, V_Handle);
   
   E_u1=M(FreeDofs,FreeDofs)\A_u1(FreeDofs,FreeDofs); E_u1=A_u1(FreeDofs,FreeDofs)'*E_u1;
   E_u2=M(FreeDofs,FreeDofs)\A_u2(FreeDofs,FreeDofs); E_u2=A_u2(FreeDofs,FreeDofs)'*E_u2;
   E_u3=M(FreeDofs,FreeDofs)\A_u3(FreeDofs,FreeDofs); E_u3=A_u3(FreeDofs,FreeDofs)'*E_u3;
   E_u4=M(FreeDofs,FreeDofs)\A_u4(FreeDofs,FreeDofs); E_u4=A_u4(FreeDofs,FreeDofs)'*E_u4;
   E_u5=M(FreeDofs,FreeDofs)\A_u5(FreeDofs,FreeDofs); E_u5=A_u5(FreeDofs,FreeDofs)'*E_u5;
   E_u6=M(FreeDofs,FreeDofs)\A_u6(FreeDofs,FreeDofs); E_u6=A_u6(FreeDofs,FreeDofs)'*E_u6;
   E_s=M(FreeDofs,FreeDofs)\A_s(FreeDofs,FreeDofs); E_s=A_s(FreeDofs,FreeDofs)'*E_s;
   E_df=M(FreeDofs,FreeDofs)\A_df(FreeDofs,FreeDofs); E_df=A_df(FreeDofs,FreeDofs)'*E_df;
   E_supg=M(FreeDofs,FreeDofs)\A_supg(FreeDofs,FreeDofs); E_supg=A_supg(FreeDofs,FreeDofs)'*E_supg;
   
%    err(i,1) =min(eig(full(E_up(FreeDofs,FreeDofs)),full(B(FreeDofs,FreeDofs)))); 
%    err(i,2) =min(eig(full(E_st(FreeDofs,FreeDofs)),full(B(FreeDofs,FreeDofs))));
%    err(i,3) =min(eig(full(E_supg(FreeDofs,FreeDofs)),full(B(FreeDofs,FreeDofs))));
   
   err(i,1) = 1/eigs(inv(full(E_u1))*H(FreeDofs,FreeDofs),[],1,'lm'); 
   err(i,2) = 1/eigs(inv(full(E_u2))*H(FreeDofs,FreeDofs),[],1,'lm');
   err(i,3) = 1/eigs(inv(full(E_u3))*H(FreeDofs,FreeDofs),[],1,'lm');
   err(i,4) = 1/eigs(inv(full(E_u4))*H(FreeDofs,FreeDofs),[],1,'lm');
   %err(i,5) = eigs(E_u5(FreeDofs,FreeDofs),H(FreeDofs,FreeDofs),1,'sm');
   %err(i,6) = eigs(E_u6(FreeDofs,FreeDofs),H(FreeDofs,FreeDofs),1,'sm');

   err(i,5) = eigs(E_s,H(FreeDofs,FreeDofs),1,'sm');
   err(i,6) = eigs(E_df,H(FreeDofs,FreeDofs),1,'sm');
   %err(i,6) = eigs(E_supg,H(FreeDofs,FreeDofs),1,'sm');
   
   Dofs(i)=size(H,1);
   
 end; 
  fig = figure('Name','Discrete Inf-Sup Condition');
  plot(Dofs,sqrt(err(:,4)),'go-',Dofs,sqrt(err(:,1)),'ro-',...
      Dofs,sqrt(err(:,2)),'bo-',Dofs,sqrt(err(:,3)),'yo-',...
      Dofs,sqrt(err(:,5)),'co-',Dofs,sqrt(err(:,6)),'mo-'); grid('on');
  set(gca,'XScale','log','YScale','log');
  xlabel('{\bf DOFS}');
  ylabel('{\bf infsup-constant}');
  legend('up1','up2','up3',...
      'up4','standard','supg','Location','West');
%   p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,4)),1);
%   add_Slope(gca,'West',p(1),'g-');
%   p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,1)),1);
%   add_Slope(gca,'SouthEast',p(1),'r-');
%   p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,2)),1);
%   add_Slope(gca,'East',p(1),'b-');
%   p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,3)),1);
%   add_Slope(gca,'Southwest',p(1),'y-');
%   p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,5)),1);
%   add_Slope(gca,'NorthWest',p(1),'c-');
%   p = polyfit(log(Dofs(NREFS-3:NREFS)),log(err(NREFS-3:NREFS,6)),1);
%   add_Slope(gca,'NorthEast',p(1),'m-');
save ('data_QFE', 'err', 'Dofs')
print -depsc 'infsup_QFE.eps'
clear all;
