% Ratio between energy error and estimators on regular refined meshes. 
% The script generates the figures:
%  ErrEst_RES_REG_3.eps,  ErrEst_REC_REG_3.eps,  ErrEst_HIER_REG_3.eps,
%  ErrEst_RES_REG_4.eps,  ErrEst_REC_REG_4.eps,  ErrEst_HIER_REG_4.eps,
%  ratio_ErrEst_REG.eps.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize constants
  
  F_HANDLE = @f_LShap;           % Right hand side source term
  GD_HANDLE = @g_D_LShap;        % Dirichlet boundary data
  UEX_HANDLE = @grad_uex_LShap;  % Exact solution
  THETA = .5;                    % Refinement strategy
  NREFS = 8;                     % Number of red refinement steps
  PLOTS = [5];                 % Plots of error distributions       
   
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat'); 
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);         
  Mesh = add_Edges(Mesh);                                
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
  Mesh.BdFlags(Loc) = -1;
  
  % Compute the error and the error estimators
  
  for iter = 1:NREFS
  
    nDofs(iter) = size(Mesh.Coordinates,1); 
   
    % Compute FE solution
  
    A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
    L = assemLoad_LFE(Mesh,P7O6(),F_HANDLE);
    [U,FreeDofs] = assemDir_LFE(Mesh,-1,GD_HANDLE);
    L = L - A*U;
    if(~isempty(FreeDofs))
      U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    end
    
    % Compute RED discrete errors
    
    [Err(iter),ErrDist] = H1SErr_LFE(Mesh,U,P7O6(),UEX_HANDLE);
    ErrDist = sqrt(ErrDist);
    
    % Each estimated error
    
    Mesh = add_Edge2Elem(Mesh);
    Mesh = add_Patches(Mesh);
    
    Eta_RES_K = ErrEst_RES(U,Mesh,P7O6(),F_HANDLE);    
    Eta_RES(iter) = sqrt(sum(Eta_RES_K.^2));
    
    Eta_REC_K = ErrEst_REC(U,Mesh,P7O6());    
    Eta_REC(iter) = sqrt(sum(Eta_REC_K.^2));
    
    Eta_HIER_K = ErrEst_HIER(U,Mesh,@Res_Lapl,@STIMA_Lapl_ErrEst,P7O6(),F_HANDLE);  
    Eta_HIER(iter) = sqrt(sum(Eta_HIER_K.^2));
    
    % Generate plots of error distributions
    
    if(~isempty(find(iter == PLOTS)))
        
      % Plot distribution of true error
      fig = plot_EST(ErrDist,Mesh);
      colorbar;
      title(['{\bf Distribution of true error, refinement level ' int2str(iter) ' }'])
      xlabel(['{\bf' sprintf('error = %2.2e',Err(iter)) '}']);
      FileName = ['ErrDist_REG_' int2str(iter) '.eps'];
      print('-depsc', FileName);
      
        % Generate figure
        
        fig = plot_EST(Eta_RES_K,Mesh);
        colorbar;
        title(['{\bf Residual based estimate, refinement level ' int2str(iter) ' }']);
        xlabel(['{\bf' sprintf('Eta = %2.2e',Eta_RES(iter)) '}']);
        text(.1,-.4,'# Vertices','FontWeight','bold');
        text(.6,-.4,':','FontWeight','bold');
        text(.7,-.4,int2str(size(Mesh.Coordinates,1)),'FontWeight','bold');
        text(.1,-.55,'# Elements','FontWeight','bold');
        text(.6,-.55,':','FontWeight','bold');
        text(.7,-.55,int2str(size(Mesh.Elements,1)),'FontWeight','bold');
        text(.1,-.7,'# Edges','FontWeight','bold');    
        text(.6,-.7,':','FontWeight','bold');
        text(.7,-.7,int2str(size(Mesh.Edges,1)),'FontWeight','bold');  
      
        % Print  to .eps file
        
        FileName = ['ErrEst_RES_REG_' int2str(iter) '.eps'];
        print('-depsc', FileName);
%        close(fig);
%        system(['gv ' FileName ' &']);
                
        % Generate figure
        
        fig = plot_EST(Eta_REC_K,Mesh);
        colorbar;
        title(['{\bf Recovery based estimate, refinement level ' int2str(iter) ' }']);
        xlabel(['{\bf' sprintf('Eta = %2.2e',Eta_REC(iter)) '}']);
        text(.1,-.4,'# Vertices','FontWeight','bold');
        text(.6,-.4,':','FontWeight','bold');
        text(.7,-.4,int2str(size(Mesh.Coordinates,1)),'FontWeight','bold');
        text(.1,-.55,'# Elements','FontWeight','bold');
        text(.6,-.55,':','FontWeight','bold');
        text(.7,-.55,int2str(size(Mesh.Elements,1)),'FontWeight','bold');
        text(.1,-.7,'# Edges','FontWeight','bold');    
        text(.6,-.7,':','FontWeight','bold');
        text(.7,-.7,int2str(size(Mesh.Edges,1)),'FontWeight','bold');  
      
        % Print  to .eps file
        
        FileName = ['ErrEst_REC_REG_' int2str(iter) '.eps'];
        print('-depsc', FileName);
%        close(fig);
%        system(['gv ' FileName ' &']);
                
        % Generate figure
        
        fig = plot_EST(Eta_HIER_K,Mesh);
        colorbar;
        title(['{\bf Hierarchical estimate, refinement level ' int2str(iter) ' }']);
        xlabel(['{\bf' sprintf('Eta = %2.2e',Eta_HIER(iter)) '}']);
        text(.1,-.4,'# Vertices','FontWeight','bold');
        text(.6,-.4,':','FontWeight','bold');
        text(.7,-.4,int2str(size(Mesh.Coordinates,1)),'FontWeight','bold');
        text(.1,-.55,'# Elements','FontWeight','bold');
        text(.6,-.55,':','FontWeight','bold');
        text(.7,-.55,int2str(size(Mesh.Elements,1)),'FontWeight','bold');
        text(.1,-.7,'# Edges','FontWeight','bold');    
        text(.6,-.7,':','FontWeight','bold');
        text(.7,-.7,int2str(size(Mesh.Edges,1)),'FontWeight','bold');  
      
        % Print  to .eps file
        
        FileName = ['ErrEst_HIER_REG_' int2str(iter) '.eps'];
        print('-depsc', FileName);
 %       close(fig);
 %       system(['gv ' FileName ' &']);
        
    end
    
    % Generate new mesh
    
    Mesh = refine_REG(Mesh);  
    
    % Print out information
    
    fprintf('Err (RED)  :  %2.4e,   iter  :  %d\n',Err(iter),iter);
  
  end
  
  % Generate figures 
  
  fig = figure();
  plot(nDofs,Eta_RES./Err,'-x', ...
       nDofs,Eta_REC./Err,'-*', ...
       nDofs,Eta_HIER./Err,'->');
  grid('on');
  set(gca,'XScale','log');
  title('{\bf Ratio between energy error and estimators}');
  xlabel('{\bf Dofs [log]}');
  ylabel('{\bf Ratio}');
  legend(['Err. est. RES reg. ref.'], ...
         ['Err. est. REC reg. ref.'], ...
         ['Err. est. HIER reg. ref.'], ...
         'Location','East');
  print('-depsc','ratio_ErrEst_REG.eps');
  % close(fig);
  %!gv ratio_ErrEst_REG.eps &

  % Clear memory
  
  clear all;