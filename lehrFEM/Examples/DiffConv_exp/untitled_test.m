% Initialize constants
  clear Mesh
 
  JIG=1;
  d1=0.1;                    %impact of supg modification
  d2=0.1;
  a=10^0;                %amount of diffusivity
  NREFS =5;
  % select test code
  d=getData(7);
     
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
   %
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
   
   % old_version
   weights1=[0 0 0 1/3 1/3 1/3 0];
   weights2=[1/12 1/12 1/12 0 0 0 3/4];
   weights3=[0.05 0.05 0.05 4/30 4/30 4/30, 0.45];
   weights4=[1/3 1/3 1/3 0 0 0 0];
   [M1,UP1,C1]=assemMat_UpQFE(New_Mesh,d.V_Handle,weights1);
   [M2,UP2,C2]=assemMat_UpQFE(New_Mesh,d.V_Handle,weights2);
   [M3,UP3,C3]=assemMat_UpQFE(New_Mesh,d.V_Handle,weights3);
   [M4,UP4,C4]=assemMat_UpQFE(New_Mesh,d.V_Handle,weights4);
   A1=M1*UP1;
   A2=M2*UP2+C2;
   A3=M3*UP3+C3;
   A4=M4*UP4+C4;
   
   
   % new version
   % weight1
   [qr_bnd qr_inn]=split_QuadRule(Duffy(TProd(gaulob(0,1,10))));
   if (~isempty(qr_inn.x))
        newA1=assemMat_QFE(New_Mesh,@STIMA_Conv_QFE,d.V_Handle,qr_inn);
        if (~isempty(qr_bnd.x))
            newA1=newA1+assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   else
        if (~isempty(qr_bnd.x))
            newA1=assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   end
   
   
   % weight2
   [qr_bnd qr_inn]=split_QuadRule(P4O3());
    if (~isempty(qr_inn.x))
        newA2=assemMat_QFE(New_Mesh,@STIMA_Conv_QFE,d.V_Handle,qr_inn);
        if (~isempty(qr_bnd.x))
            newA2=newA2+assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   else
        if (~isempty(qr_bnd.x))
            newA2=assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   end
   
   % weight3
   [qr_bnd qr_inn]=split_QuadRule(P7O4());
    if (~isempty(qr_inn.x))
        newA3=assemMat_QFE(New_Mesh,@STIMA_Conv_QFE,d.V_Handle,qr_inn);
        if (~isempty(qr_bnd.x))
            newA3=newA3+assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   else
        if (~isempty(qr_bnd.x))
            newA3=assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   end
   
   % weight4
   [qr_bnd qr_inn]=split_QuadRule(P3O2());
   if (~isempty(qr_inn.x))
        newA4=assemMat_QFE(New_Mesh,@STIMA_Conv_QFE,d.V_Handle,qr_inn);
        if (~isempty(qr_bnd.x))
            newA4=newA4+assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   else
        if (~isempty(qr_bnd.x))
            newA4=assemMat_UPQFE2(New_Mesh,@STIMA_UPQFE2,d.V_Handle,qr_bnd);
        end
   end
     
  err(i,1) = norm(full(A1-newA1));
  err(i,2) = norm(full(A2-newA2));
  err(i,3) = norm(full(A3-newA3));
  err(i,4) = norm(full(A4-newA4));
  Dofs(i)=size(New_Mesh.Coordinates,1); 
  Mesh = refine_REG(Mesh);
  
 end; 
  fig = figure('Name','Discretization error');
  plot(h,err(:,4),'go-',h,err(:,1),'ro-',h,err(:,2),'bo-',h,err(:,3),'yo-'); grid('on');
  set(gca,'XScale','log','YScale','log');
  xlabel('{\bf h}');
  ylabel('{\bf Error}');
  legend('L^2-error u (up1)','L^2-error u (up2)','L^2-error u (up3)',...
      'L^2-error u (up4)','Location','NorthEast');
clear all;