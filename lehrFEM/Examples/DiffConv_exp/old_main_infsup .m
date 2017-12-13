% driver routine for several diffusion convection discretizations
% see below
   clear Mesh
  % Initialize constants
  JIG=1;
  d1=1;                      %impact of supg modification
  d2=1;
  a=10^-6               %amount of diffusivity
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
 % Mesh = refine_REG(Mesh);
 
 tic;
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
   
   M = assemMat_LFE(NewMesh,@MASS_LFE);
  
   %diagonal Mass Matrix
   
   MassZero = assemMat_Mass0fD(NewMesh); 
   
   %Convection term 
   
   % upwinding Quadrature Rule= discrete differential forms
   
   D_up = assemMat_LFE(NewMesh, @STIMA_ContrGrad_Up, V_Handle);
   D_up = MassZero * D_up;
   % or equivalently
                  TopGrad=assemMat_TopGrad(NewMesh);
   %               ContrOne=assemMat_Contr1f(NewMesh,d.V_Handle);
   %               D=ContrOne*TopGrad;           
      
   % standard Galerkin
   
   D_st = assemMat_LFE(NewMesh,@STIMA_Conv_LFE ,V_Handle,QuadRule);
   
   % SUPG
   
   D_supg2 = assemMat_LFE(NewMesh,@STIMA_SUPG_LFE,QuadRule, V_Handle,a,d1,d2);
   D_supg  = D_st+D_supg2;
   
   %selfcooked
   M_QFELFE=assemMat_MASS_QFELFE(NewMesh,P3O3());
   ContrOneQFE=assemMat_Contr1f_QFE(NewMesh,V_Handle);  % contraction of one forms
   D_c=M_QFELFE*ContrOneQFE*TopGrad;
   
   %Laplace + convection
   
   A_up = (a * A + D_up);
   A_st = (a* A + D_st);
   A_supg = (a * A + D_supg);
   A_c = (a * A + D_c);
      
   % discrete inf-sup-condition
   
   FreeDofs = FreeDofs_LFE(NewMesh,[-1,-2]);
   
   B=assemMat_LFE_BdLayer(NewMesh,@STIMA_InfSup_LFE,QuadRule, V_Handle,1,0,0);
   E_up = M(FreeDofs,FreeDofs)\A_up(FreeDofs,FreeDofs); E_up=A_up(FreeDofs,FreeDofs)'*E_up;
   E_st = M(FreeDofs,FreeDofs)\A_st(FreeDofs,FreeDofs); E_st=A_st(FreeDofs,FreeDofs)'*E_st;
   E_supg=M(FreeDofs,FreeDofs)\A_supg(FreeDofs,FreeDofs); E_supg=A_supg(FreeDofs,FreeDofs)'*E_supg;
   E_c = M(FreeDofs,FreeDofs)\A_c(FreeDofs,FreeDofs); E_c=A_c(FreeDofs,FreeDofs)'*E_c;   
%    err(i,1) =min(eig(full(E_up(FreeDofs,FreeDofs)),full(B(FreeDofs,FreeDofs)))); 
%    err(i,2) =min(eig(full(E_st(FreeDofs,FreeDofs)),full(B(FreeDofs,FreeDofs))));
%    err(i,3) =min(eig(full(E_supg(FreeDofs,FreeDofs)),full(B(FreeDofs,FreeDofs))));
   
%    
   err(i,1) = 1/eigs(inv(full(E_up))*B(FreeDofs,FreeDofs),[],1,'LM'); 
   err(i,2) = 1/eigs(inv(full(E_st))*B(FreeDofs,FreeDofs),[],1,'LM');
   err(i,3) = 1/eigs(inv(full(E_supg))*B(FreeDofs,FreeDofs),[],1,'LM');
   err(i,4) = 1/eigs(inv(full(E_c))*B(FreeDofs,FreeDofs),[],1,'LM');
    
%    err(i,1) = power_iter_infsup(A_up(FreeDofs,FreeDofs), M(FreeDofs,FreeDofs), B(FreeDofs,FreeDofs), 10^-9,10^5)
%    err(i,2) = power_iter_infsup(A_st(FreeDofs,FreeDofs), M(FreeDofs,FreeDofs), B(FreeDofs,FreeDofs), 10^-9,10^5)
%    err(i,3) = power_iter_infsup(A_supg(FreeDofs,FreeDofs), M(FreeDofs,FreeDofs), B(FreeDofs,FreeDofs), 10^-9,10^5)
%    err(i,4) = power_iter_infsup(A_c(FreeDofs,FreeDofs), M(FreeDofs,FreeDofs), B(FreeDofs,FreeDofs), 10^-9,10^5)
   
   
   Dofs(i)=size(NewMesh.Coordinates,1);
   
   % plot_LFE(U_up,NewMesh); colorbar;
   % plot_LFE(U_st,NewMesh);  colorbar; 
   % plot_LFE(U_supg,NewMesh); colorbar;

   % plot_LFE(U_c,NewMesh); colorbar;
  end; 
  toc
  
  fig = figure('Name','Discret Inf-Sup-Condition');
  plot(Dofs,sqrt(err(:,1)),'r--',...
	Dofs,sqrt(err(:,2)),'b-.',...
	Dofs,sqrt(err(:,3)),'g-',...
	Dofs,sqrt(err(:,4)),'c-'); grid('on');
  set(gca,'XScale','log','YScale','log');
  xlabel('{\bf Dofs}');
  ylabel('{\bf Inf-Sup-Constant}');
  legend('Up','St','SUPG','Location','NorthEast')
  save ('data_LFE', 'err', 'Dofs')
  print -depsc 'infsup_LFE.eps'
clear all;
