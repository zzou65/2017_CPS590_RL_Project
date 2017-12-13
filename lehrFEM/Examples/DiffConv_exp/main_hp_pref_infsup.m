% Runs script for convergence rates of hp-FEM, p refinement

% Copyright 2006-2006 Patrick Meury
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize constants
  clear Mesh;
  clear Mesh_hp;
  NREFS =7;                           % Number of p refinements  
  a=10^-6; 
  P_Add_Deg=3;
  
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
  
  CNodes = transpose(1:size(Mesh.Coordinates,1));
  
  % Prepare mesh for longest edge bisection 
  
  nDofs = zeros(1,NREFS);
  H1S_error = zeros(1,NREFS);
  L2_error = zeros(1,NREFS);
 
  Mesh = init_LEB(Mesh);
 % Mesh = refine_REG(Mesh);
  for i = 2:NREFS
    P_DEG=i;
    
    %Mesh = refine_hp(Mesh,CNodes);
  %  Mesh = refine_REG(Mesh);
    
    % Generate mesh data structure for hp-FEM
  
    Mesh_hp.Coordinates = Mesh.Coordinates;
    Mesh_hp.Elements = Mesh.Elements;
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
    
    %[EDofs,CDofs] = assign_pdeg(Mesh_hp,CNodes,P_DEG); 
    %Elem2Dof = build_DofMaps(Mesh_hp,EDofs,CDofs);
    % pmax = max(CDofs)+1;
    
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
    %velocity=V_Handle(QuadRule_2D.x);
    
    %Global stiffness matrix
    M = assemMat_hp(Mesh_hp,Elem2Dof,@MASS_hp,QuadRule_2D,Shap_2D);
    
    A = assemMat_hp(Mesh_hp,Elem2Dof,@STIMA_Lapl_hp,QuadRule_2D,Shap_2D); 
    
    B = assemMat_hp(Mesh_hp,Elem2Dof,@STIMA_Conv_hp,QuadRule_2D,Shap_2D,V_Handle);
    
    [U,FreeDofs] = assemDir_hp(Mesh_hp,Elem2Dof,[-1, -2],QuadRule_1D,Shap_1D,@(x,varargin)1);
    
    H=assemMat_hp_BdLayer(Mesh_hp,Elem2Dof,@STIMA_InfSup_HP,QuadRule_2D, Shap_2D, V_Handle);
    
    h(i)=get_MeshWidth(Mesh);
    nDofs(i) = size(U,1);
  
    % new version
    for p=1:P_Add_Deg
    
        [qr_bnd qr_inn]=split_QuadRule(Duffy(TProd(gaulob(0,1,2*pmax+p-1))));
        %[qr_bnd qr_inn]=split_QuadRule(stab_QuadRule);
        if (~isempty(qr_inn.x))
            Shap_inn = shap_hp(qr_inn.x,pmax);
            B_stab=assemMat_hp(Mesh_hp,Elem2Dof,@STIMA_Conv_hp,qr_inn,Shap_inn,V_Handle);
            if (~isempty(qr_bnd.x))
                B_stab=B_stab+assemMat_Uphp(Mesh_hp, Elem2Dof, V_Handle, qr_bnd, Shap_2D,pmax);   
            end
        else
            if (~isempty(qr_bnd.x))
                B_stab=assemMat_Uphp(Mesh_hp, Elem2Dof, V_Handle, qr_bnd, Shap_2D,pmax);   
            end
        end

        S_stab=a*A+B_stab;
        
        % discrete inf-sup-condition
    
        E_stab=M\S_stab; E_stab=S_stab'*E_stab;
   
        err(p,i) = eigs(E_stab(FreeDofs,FreeDofs),H(FreeDofs,FreeDofs),1,'sm'); 
        
    end

    S=a*A+B;
    E=M\S; E=S'*E;
    err(p+1,i) = eigs(E(FreeDofs,FreeDofs),H(FreeDofs,FreeDofs),1,'sm');
  end
    
  % Generate figure
  
  fig = figure();
  plot(1:NREFS,err); hold on;
  title('{\bf Discretization errors for hp-FEM}')
  xlabel('{\bf # pol. degree}');
  ylabel('{\bf infsup}');
  set(gca,'YScale','log');
  legend('stab','standard ','Location','NorthEast')
  print -depsc 'hp_p_infsup.eps'
  % Clear memory
  
  clear all;
  