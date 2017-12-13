function main_hp_href_infsup(degree)
% Runs script for convergence rates of hp-FEM. h-refinement

% Copyright 2006-2006 Patrick Meury
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

  % Initialize constants
  JIG=2;
  clear Mesh;
  clear Mesh_hp;
  P_DEG = degree;                            % polynomial approximation
  NREFS = 4;                           % Number of mesh refinements  
  a = 10^-8;
  P_Add_Deg=1;
  %stab_QuadRule=Duffy(TProd(gaulob(0,1,2*P_DEG)))
  %stab_QuadRule=P10O5();
  V_Handle=...
            @(x,varargin)[ones(size(x,1),1) ones(size(x,1),1)];  
 
  QuadRule = P7O6(); 
  Mesh.Coordinates =[-1 -1; 1 -1; 1 1; -1 1]; 
  Mesh.Elements = [1 2 3;1 3 4];
  Mesh = add_Edges(Mesh);
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
  Mesh.BdFlags(Loc) = [-1 -1 -2 -2];   % -1 Inflow
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
  
  CNodes = transpose(1:size(Mesh.Coordinates,1));
  
  % Prepare mesh for longest edge bisection 
  
  nDofs = zeros(1,NREFS);
  h = zeros(1,NREFS);
  err = zeros(P_Add_Deg+1,NREFS);
  
  Mesh = init_LEB(Mesh);
  Mesh = refine_REG(Mesh);
  %for i = 1:NREFS
  i=1; dofs=0;
  while (i<=NREFS)
      i
%     %Mesh = refine_hp(Mesh,CNodes);
%     if (P_DEG==1) 
%          Mesh = refine_REG(Mesh);
%     end
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
  
    Mesh_hp.Coordinates = Mesh.Coordinates;
    Mesh_hp.Elements = Mesh.Elements;
    Mesh_hp.ElemFlag = zeros(size(Mesh_hp.Elements,1),1);
    Mesh_hp = add_Edges(Mesh_hp);
    if(isfield(Mesh_hp,'BdFlags'))
      Mesh_hp = rmfield(Mesh_hp,'BdFlags');
    end
    Loc = get_BdEdges(Mesh_hp);
    Mesh_hp.BdFlags = zeros(size(Mesh_hp.Edges,1),1);
    Mesh_hp.BdFlags(Loc) = Mesh.BdFlags(Loc);
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
    
    %Global stiffness matrix
    M = assemMat_hp(Mesh_hp,Elem2Dof,@MASS_hp,QuadRule_2D,Shap_2D);
    
    A = assemMat_hp(Mesh_hp,Elem2Dof,@STIMA_Lapl_hp,QuadRule_2D,Shap_2D); 
       
    [U,FreeDofs] = assemDir_hp(Mesh_hp,Elem2Dof,[-1, -2],QuadRule_1D,Shap_1D,@(x,varargin)0);
    
    H=assemMat_hp_BdLayer(Mesh_hp,Elem2Dof,@STIMA_InfSup_HP,QuadRule_2D, Shap_2D, V_Handle);
    
    nDofs(i) = size(U,1);
    dofs=nDofs(i)
  
    % new version
    for p=1:P_Add_Deg
    
        [qr_bnd qr_inn]=split_QuadRule(Duffy(TProd(gaulob(0,1,pmax+p))));
        %[qr_bnd qr_inn]=split_QuadRule(P7O4());
        if (~isempty(qr_inn.x))
            Shap_inn = shap_hp(qr_inn.x,pmax);
            B=assemMat_hp(Mesh_hp,Elem2Dof,@STIMA_Conv_hp,qr_inn,Shap_inn,V_Handle);
            if (~isempty(qr_bnd.x))
                Shap_bnd = shap_hp(qr_bnd.x,pmax);
                B=B+assemMat_Uphp(Mesh_hp, Elem2Dof, V_Handle, qr_bnd, Shap_bnd,pmax);   
            end
        else
            if (~isempty(qr_bnd.x))
                Shap_bnd = shap_hp(qr_bnd.x,pmax);
                B=assemMat_Uphp(Mesh_hp, Elem2Dof, V_Handle, qr_bnd, Shap_bnd,pmax);
            end
        end

        S=a*A+B;
        
        % discrete inf-sup-condition
    
        E=M(FreeDofs,FreeDofs)\S(FreeDofs,FreeDofs); E=S(FreeDofs,FreeDofs)'*E;
        clear U S B;
        err(p,i) = 1/eigs(inv(full(E))*H(FreeDofs,FreeDofs),[],1,'lm');
        
    end

    B = assemMat_hp(Mesh_hp,Elem2Dof,@STIMA_Conv_hp,QuadRule_2D,Shap_2D,V_Handle);
 
    S=a*A+B;
    E=M(FreeDofs,FreeDofs)\S(FreeDofs,FreeDofs); E=S(FreeDofs,FreeDofs)'*E;
    clear S B M ;
    err(p+1,i) = 1/eigs(inv(full(E))*H(FreeDofs,FreeDofs),[],1,'lm');
    
    Mesh = refine_REG(Mesh);
    i=i+1;
  end
    
  % Generate figure
  
  fig = figure();
  plot(nDofs,err,'*-'); hold on;
  xlabel('{\bf # DOFs}');
  ylabel('{\bf infsup}');
  set(gca,'XScale','log','YScale','log');
  l_matrix=int2str((1:P_Add_Deg+1)');
  
  legend(l_matrix,'Location','East')
  filename=['hp_h_infsup_p',int2str(degree),'.eps'];
  print('-depsc' , filename);
  % Clear memory
  
  clear all;