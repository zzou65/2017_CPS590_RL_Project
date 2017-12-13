% Runs script for convergence rates of hp-FEM. h-refinement

% Copyright 2006-2006 Patrick Meury
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize constants
  JIG=1;
  clear Mesh;
  clear Mesh_hp;
  P_DEG =3;                            % polynomial approximation
  NREFS =6;                           % Number of mesh refinements  
  a=10^0;  
  d=getData(1);
  
  % Initialize mesh
   
  Mesh.Coordinates =d.Coordinates; 
  Mesh.Elements = d.Elements;
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = d.boundtype;
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);  
  
  CNodes = transpose(1:size(Mesh.Coordinates,1));
  
  % Prepare mesh for longest edge bisection 
  
  nDofs = zeros(1,NREFS);
  h = zeros(1,NREFS);
  H1S_error = zeros(1,NREFS);
  L2_error = zeros(1,NREFS);
 
  Mesh = init_LEB(Mesh);
  % Mesh = refine_REG(Mesh);
  for i = 1:NREFS
    
     %Mesh = refine_hp(Mesh,CNodes);
     Mesh = refine_REG(Mesh);
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
    
    % Generate mesh data structure for hp-FEM
  
    Mesh_hp.Coordinates = New_Mesh.Coordinates;
    Mesh_hp.Elements = New_Mesh.Elements;
    Mesh_hp.ElemFlag = zeros(size(Mesh_hp.Elements,1),1);
    Mesh_hp = add_Edges(Mesh_hp);
    if(isfield(Mesh_hp,'BdFlags'))
      Mesh_hp = rmfield(Mesh_hp,'BdFlags');
    end
    Loc = get_BdEdges(Mesh_hp);
    Mesh_hp.BdFlags = zeros(size(Mesh_hp.Edges,1),1);
    Mesh_hp.BdFlags(Loc) = -1;
    Mesh_hp = add_Edge2Elem(Mesh_hp);
  
    % Assign polynomial degrees and build dof maps
   
    % inhomogen
    
    %     [EDofs,CDofs] = assign_pdeg(Mesh_hp,CNodes,NREFS); 
    %     Elem2Dof = build_DofMaps(Mesh_hp,EDofs,CDofs);
    %     pmax = max(CDofs)+1;
    
    % homogen
    
    EDofs = (P_DEG-1)*ones(size(Mesh_hp.Edges,1),1);
    if(P_DEG > 2)
      CDofs = (P_DEG-1)*(P_DEG-2)/2*ones(size(Mesh_hp.Elements,1),1);
    else
      CDofs = zeros(size(Mesh_hp.Elements,1),1);
    end
    Elem2Dof = build_DofMaps(Mesh_hp,EDofs,CDofs);  
    pmax = P_DEG;
    
    % Build shape functions and quadrature rules

    QuadRule_1D = gauleg(0,1,2*pmax);
    Shap_1D = shap_hp([QuadRule_1D.x zeros(size(QuadRule_1D.x))],pmax);
    QuadRule_2D = Duffy(TProd(QuadRule_1D));
    Shap_2D = shap_hp(QuadRule_2D.x,pmax);

    QuadRule_1D_quad = gauleg(0,1,10*pmax);
    Shap_1D_quad = shap_hp([QuadRule_1D_quad.x zeros(size(QuadRule_1D_quad.x))],pmax);
    QuadRule_2D_quad = Duffy(TProd(QuadRule_1D_quad));
    Shap_2D_quad = shap_hp(QuadRule_2D_quad.x,pmax);
    
    %Global stiffness matrix
    
    A = assemMat_hp(Mesh_hp,Elem2Dof,@STIMA_Lapl_hp,QuadRule_2D,Shap_2D); 
    
    B = assemMat_hp(Mesh_hp,Elem2Dof,@STIMA_Conv_hp,QuadRule_2D,Shap_2D,d.V_Handle);
    
    % new version
    %[qr_bnd qr_inn]=split_QuadRule(Duffy(TProd(gaulob(0,1,2*pmax+1))));
    [qr_bnd qr_inn]=split_QuadRule(NCC(2*pmax));
    %[qr_bnd qr_inn]=split_QuadRule(P10O4());
    %[qr_bnd qr_inn]=split_QuadRule(Duffy(TProd(gauleg(0,1,3))));
    if (~isempty(qr_inn.x))
        Shap_inn = shap_hp(qr_inn.x,pmax);
        B_stab=assemMat_hp(Mesh_hp,Elem2Dof,@STIMA_Conv_hp,qr_inn,Shap_inn,d.V_Handle);
        if (~isempty(qr_bnd.x))
            Shap_bnd = shap_hp(qr_bnd.x,pmax);
            B_stab=B_stab+assemMat_Uphp(Mesh_hp, Elem2Dof, d.V_Handle, qr_bnd, Shap_bnd,pmax);   
        end
    else
        if (~isempty(qr_bnd.x))
            Shap_bnd = shap_hp(qr_bnd.x,pmax);
            B_stab=assemMat_Uphp(Mesh_hp, Elem2Dof, d.V_Handle, qr_bnd, Shap_bnd,pmax);   
        end
    end

    L = assemLoad_hp(Mesh_hp,Elem2Dof,QuadRule_2D_quad,Shap_2D_quad,d.SOL_Handle,a);
  
    % Incoporate Dirichlet boundary conditions
    
    [U,FreeDofs] = assemDir_hp(Mesh_hp,Elem2Dof,-1,QuadRule_1D_quad,Shap_1D_quad,d.U_EX_Handle,a);
    U_stab=U;
    S=a*A+B;
    S_stab=a*A+B_stab;
    L_stab = L - S_stab*U;
    L = L - S*U;
    
    % Solve the linear system
  
    U_stab(FreeDofs) = S_stab(FreeDofs,FreeDofs)\L_stab(FreeDofs);
    U(FreeDofs) = S(FreeDofs,FreeDofs)\L(FreeDofs);
    
    % Compute discretization errors
    
    h(i)=get_MeshWidth(Mesh);
    nDofs(i) = size(U,1);
    H1S_error(i) = H1SErr_hp(Mesh_hp,U_stab,Elem2Dof,QuadRule_2D_quad,Shap_2D_quad,d.GRAD_U_EX_Handle,0,a);
    L2_error(i) = H1SErr_hp(Mesh_hp,U,Elem2Dof,QuadRule_2D_quad,Shap_2D_quad,d.GRAD_U_EX_Handle,0,a);
    H1S_error2(i) = L2Err_hp(Mesh_hp,U_stab,Elem2Dof,QuadRule_2D,Shap_2D,d.U_EX_Handle,0,a);
    L2_error2(i) = L2Err_hp(Mesh_hp,U,Elem2Dof,QuadRule_2D,Shap_2D,d.U_EX_Handle,0,a);
    
   % plot_hp(U,Mesh_hp,Elem2Dof,pmax);
   % plot_hp_vtk(U_stab-U,Mesh_hp,Elem2Dof,pmax,'vtktest');
   % colorbar;
   % plot_hp((U_stab-U),Mesh_hp,Elem2Dof,pmax);
   % colorbar;
%     
%     fig = figure('Name','Polynomial degrees');  
%     patch('Faces',Mesh_hp.Elements, ...
%         'Vertices',Mesh_hp.Coordinates, ...
%         'FaceVertexCData',CDofs, ...
%         'EdgeColor','k', ...
%         'FaceColor','flat');
%     set(gca,'CLim',[1 NREFS],'DataAspectRatio',[1 1 1]);
%     colormap(jet);
%     alpha(.9);
%     colorbar;
%     set(gcf,'renderer','openGL');
i
  end
    
  % Generate figure
  %h=nDofs;
  
  fig = figure('Name','Convergence rates for hp-FEM');
  plot(h,H1S_error,'r-o',h,L2_error,'b-o');
  title('{\bf Discretization errors for hp-FEM}')
  xlabel('{\bf # Dofs [log]}');
  ylabel('{\bf Discretization error [log]}');
  set(gca,'XScale','log','YScale','log');
  legend('L^2-error (stab)','L^2-error ','Location','NorthEast')
  p = polyfit(log(h(NREFS-3:NREFS)),log(H1S_error(NREFS-3:NREFS)),1);
  add_Slope(gca,'SouthEast',p(1),'r-');
  p = polyfit(log(h(NREFS-3:NREFS)),log(L2_error(NREFS-3:NREFS)),1);
  add_Slope(gca,'East',p(1),'b-');
  
  fig = figure('Name','Convergence rates for hp-FEM');
  plot(h,H1S_error2,'r-o',h,L2_error2,'b-o');
  title('{\bf Discretization errors for hp-FEM}')
  xlabel('{\bf # Dofs [log]}');
  ylabel('{\bf Discretization error [log]}');
  set(gca,'XScale','log','YScale','log');
  legend('L^2-error (stab)','L^2-error ','Location','NorthEast')
  p = polyfit(log(h(NREFS-3:NREFS)),log(H1S_error2(NREFS-3:NREFS)),1);
  add_Slope(gca,'SouthEast',p(1),'r-');
  p = polyfit(log(h(NREFS-3:NREFS)),log(L2_error2(NREFS-3:NREFS)),1);
  add_Slope(gca,'East',p(1),'b-');
  print -depsc 'hp_ref.eps'
  % Clear memory
  
  clear all;