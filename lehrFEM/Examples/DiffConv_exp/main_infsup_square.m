% driver routine for several diffusion convection discretizations
% see below
   clear Mesh
  % Initialize constants
  JIG=1;
  d1=1;                      %impact of supg modification
  d2=1;
  as=10.^-[1:6]            %amount of diffusivity
  na=size(as,2);
  NREFS =3;
  % select test code
  V_Handle = ...
            @(x,varargin)[ones(size(x,1),1) ones(size(x,1),1)];
        
  QuadRule = P7O6(); 
  Mesh.Coordinates =[-1 -1; 1 -1; 1 1; -1 1]; 
  Mesh.Elements = [1 2 4;2 3 4];
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = [-1 -1 -2 -2];   % -1 Inflow
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  
  err = zeros(na, 4);
  Dofs=zeros(na, 1);
  
  Mesh = refine_REG(Mesh);
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
  end
      
      for ia=1:na
          
          a=as(ia);
          
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
          %               TopGrad=assemMat_TopGrad(NewMesh);
          %               ContrOne=assemMat_Contr1f(NewMesh,d.V_Handle);
          %               D=ContrOne*TopGrad;           
      
          % standard Galerkin
   
          D_st = assemMat_LFE(NewMesh,@STIMA_Conv_LFE ,V_Handle,QuadRule);
   
          % SUPG
   
          D_supg2 = assemMat_LFE(NewMesh,@STIMA_SUPG_LFE,QuadRule, V_Handle,a,d1,d2);
          D_supg  = D_st+D_supg2;
   
          %Laplace + convection
   
          A_up = (a * A + D_up);
          A_st = (a* A + D_st);
          A_supg = (a * A + D_supg);
      
          % discrete inf-sup-condition
   
          FreeDofs = FreeDofs_LFE(NewMesh,[-1,-2]);
   
          B=assemMat_LFE_BdLayer(NewMesh,@STIMA_InfSup_LFE,QuadRule, V_Handle,1,0,0);
          E_up = M\A_up; E_up=A_up'*E_up;
          E_st = M\A_st; E_st=A_st'*E_st;
          E_supg=M\A_supg; E_supg=A_supg'*E_supg;
   
          %    err(i,1) =min(eig(full(E_up(FreeDofs,FreeDofs)),full(B(FreeDofs,FreeDofs)))); 
          %    err(i,2) =min(eig(full(E_st(FreeDofs,FreeDofs)),full(B(FreeDofs,FreeDofs))));
          %    err(i,3) =min(eig(full(E_supg(FreeDofs,FreeDofs)),full(B(FreeDofs,FreeDofs))));
   
         err(ia,1) = eigs(E_up(FreeDofs,FreeDofs),B(FreeDofs,FreeDofs),1,'SM'); 
         err(ia,2) = eigs(E_st(FreeDofs,FreeDofs),B(FreeDofs,FreeDofs),1,'SM');
         err(ia,3) = eigs(E_supg(FreeDofs,FreeDofs),B(FreeDofs,FreeDofs),1,'SM');
   
      end
  fig = figure('Name','Discret Inf-Sup-Condition');
  plot(as,sqrt(err(:,1)),'r--',...
	as,sqrt(err(:,2)),'b-.',...
	as,sqrt(err(:,3)),'g-'); grid('on');
  set(gca,'XScale','log','YScale','log');
  xlabel('{\bf eps}');
  ylabel('{\bf Inf-Sup-Constant}');
  legend('Up','St','SUPG','Location','NorthEast')
  save ('data_LFE', 'err', 'Dofs')
  print -depsc 'infsup_LFE_diff.eps'
  
clear all;
