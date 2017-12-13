% Convergence rates and estimated errors on adaptive and regular refined
% meshes. The script generates the figures:
%  ErrEst_RES_1.eps,   ErrEst_REC_1.eps,   ErrEst_HIER_1.eps,  
%  ErrEst_RES_2.eps,   ErrEst_REC_2.eps,   ErrEst_HIER_2.eps,
%  ErrEst_RES_10.eps,  ErrEst_REC_10.eps,  ErrEst_HIER_10.eps,
%  ErrEst_RES_20.eps,  ErrEst_REC_20.eps,  ErrEst_HIER_20.eps,
%  ErrEst_RES_30.eps,  ErrEst_REC_30.eps,  ErrEst_HIER_30.eps,
%  ErrEst_RES_35.eps,  ErrEst_REC_35.eps,  ErrEst_HIER_35.eps,
%  ErrEst_RES_38.eps,  ErrEst_REC_38.eps,  ErrEst_HIER_38.eps,
%  rate_ErrEst.eps,    ratio_ErrEst.eps.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize constants
  
  F_HANDLE = @f_LShap;           % Right hand side source term
  GD_HANDLE = @g_D_LShap;        % Dirichlet boundary data
  UEX_HANDLE = @grad_uex_LShap;  % Exact solution
  THETA = .5;                    % Refinement strategy
  NITER = 38;                    % Total number of iterations
  NREFS = 8;                     % Number of red refinement steps
  PLOTS = [1 2 10 20 30 35 38];  % Plots of error distributions       
  
  % Save working environment
  
  save workspace.mat;  % Save initial constants 
  
  %-------------------------------------%
  %   Residual based error estimator    %
  %-------------------------------------%
  
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat'); 
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);         
  Mesh = add_Edges(Mesh);                                
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
  Mesh.BdFlags(Loc) = -1;
  Mesh = init_LEB(Mesh);
  
  for iter = 1:NITER
  
    nDofs(iter) = size(Mesh.Coordinates,1);  
      
    % Compute FE solution
  
    A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
    L = assemLoad_LFE(Mesh,P7O6(),F_HANDLE);
    [U,FreeDofs] = assemDir_LFE(Mesh,-1,GD_HANDLE);
    L = L - A*U;
    if(~isempty(FreeDofs))
      U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    end
    
    % Real and estimated errors
    
    Err(iter) = H1SErr_LFE(Mesh,U,P7O6(),UEX_HANDLE);
    Mesh = add_Edge2Elem(Mesh);
    Eta_K = ErrEst_RES(U,Mesh,P7O6(),F_HANDLE);    
    Eta(iter) = sqrt(sum(Eta_K.^2));
    
    % Generate plots of error distributions
    
    if(~isempty(find(iter == PLOTS)))
        
        % Generate figure
        
        fig = plot_EST(Eta_K,Mesh);
        colorbar;
        title(['{\bf Error distribution for iteration ' int2str(iter) ' (RES)}']);
        xlabel(['{\bf' sprintf('Eta = %2.2e',Eta(iter)) '}']);
        text(.1,-.4,'# Coordinates','FontWeight','bold');
        text(.6,-.4,':','FontWeight','bold');
        text(.7,-.4,int2str(size(Mesh.Coordinates,1)),'FontWeight','bold');
        text(.1,-.55,'# Elements','FontWeight','bold');
        text(.6,-.55,':','FontWeight','bold');
        text(.7,-.55,int2str(size(Mesh.Elements,1)),'FontWeight','bold');
        text(.1,-.7,'# Edges','FontWeight','bold');    
        text(.6,-.7,':','FontWeight','bold');
        text(.7,-.7,int2str(size(Mesh.Edges,1)),'FontWeight','bold');
        
        % Print  to .eps file
        
        FileName = ['ErrEst_RES_' int2str(iter) '.eps'];
        print('-depsc', FileName);
        close(fig);
        system(['gv ' FileName ' &']);
        
    end
    
    % Mark elements for refinement
    
    Eta_max = max(Eta_K);
    MarkedElem = find(Eta_K > Eta_max*THETA);
    
    % Refine mesh by largest edge bisection
    
    tmp_Mesh = Mesh;
    Mesh = refine_LEB(Mesh,MarkedElem);
    
    % Update mesh data structure
    
    Mesh = add_Edges(Mesh);
    Mesh = rmfield(Mesh,'BdFlags');
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
    Mesh.BdFlags(Loc) = -1; 
        
    % Print out information
    
    fprintf('Eta(RES)  :  %2.4e,   iter  :  %d\n',Eta(iter),iter); 
     
  end
  
  % Save data
  
  Data.nDofs(:,1) = nDofs;
  Data.Eta(:,1) = Eta;
  Data.Err(:,1) = Err;
  Data.LegendName{1} = 'RES';
  save data.mat Data;
  
  % Clear memory
  
  clear all;

  %-------------------------------------%
  %   Recovery based error estimator    %
  %-------------------------------------%
  
  % Initialize constants
  
  load workspace.mat;

  % Initialize mesh
  
  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat'); 
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);         
  Mesh = add_Edges(Mesh);                                
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
  Mesh.BdFlags(Loc) = -1;
  Mesh = init_LEB(Mesh);
      
  for iter = 1:NITER
      
    nDofs(iter) = size(Mesh.Coordinates,1);  
      
    % Compute FE solution
  
    A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
    L = assemLoad_LFE(Mesh,P7O6(),F_HANDLE);
    [U,FreeDofs] = assemDir_LFE(Mesh,-1,GD_HANDLE);
    L = L - A*U;
    if(~isempty(FreeDofs))
      U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    end
    
    % Real and estimated error
    
    Err(iter) = H1SErr_LFE(Mesh,U,P7O6(),UEX_HANDLE);  
    Mesh = add_Patches(Mesh);
    Eta_K = ErrEst_REC(U,Mesh,P7O6());    
    Eta(iter) = sqrt(sum(Eta_K.^2));
    
    % Generate plots of error distributions
    
    if(~isempty(find(iter == PLOTS)))
        
        % Generate figure
        
        fig = plot_EST(Eta_K,Mesh);
        colorbar;
        title(['{\bf Error distribution for iteration ' int2str(iter) ' (REC)}']);
        xlabel(['{\bf' sprintf('Eta = %2.2e',Eta(iter)) '}']);
        text(.1,-.4,'# Coordinates','FontWeight','bold');
        text(.6,-.4,':','FontWeight','bold');
        text(.7,-.4,int2str(size(Mesh.Coordinates,1)),'FontWeight','bold');
        text(.1,-.55,'# Elements','FontWeight','bold');
        text(.6,-.55,':','FontWeight','bold');
        text(.7,-.55,int2str(size(Mesh.Elements,1)),'FontWeight','bold');
        text(.1,-.7,'# Edges','FontWeight','bold');    
        text(.6,-.7,':','FontWeight','bold');
        text(.7,-.7,int2str(size(Mesh.Edges,1)),'FontWeight','bold');  
      
        % Print  to .eps file
        
        FileName = ['ErrEst_REC_' int2str(iter) '.eps'];
        print('-depsc', FileName);
        close(fig);
        system(['gv ' FileName ' &']);
        
    end
    
    % Mark elements for refinement
    
    Eta_max = max(Eta_K);
    MarkedElem = find(Eta_K > Eta_max*THETA);
    
    % Refine mesh by largest edge bisection
    
    tmp_Mesh = Mesh;
    Mesh = refine_LEB(Mesh,MarkedElem);
    
    % Update mesh data structure
    
    Mesh = add_Edges(Mesh);
    Mesh = rmfield(Mesh,'BdFlags');
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
    Mesh.BdFlags(Loc) = -1;   
  
    % Print out information
    
    fprintf('Eta(REC)  :  %2.4e,   iter  :  %d\n',Eta(iter),iter);
    
  end
  
  % Update data
 
  load data.mat;
  Data.nDofs(:,2) = nDofs;
  Data.Eta(:,2) = Eta;
  Data.Err(:,2) = Err;
  Data.LegendName{2} = 'REC';
  save data.mat Data;
  
  % Clear memory
  
  clear all;

  %-------------------------------------%
  %   Hierarchical error estimator      %
  %-------------------------------------%
  
  % Initialize constants
  
  load workspace.mat;
 
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat'); 
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);         
  Mesh = add_Edges(Mesh);                                
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
  Mesh.BdFlags(Loc) = -1;
  Mesh = init_LEB(Mesh);
      
  for iter = 1:NITER
  
    nDofs(iter) = size(Mesh.Coordinates,1);  
      
    % Compute FE solution
  
    A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
    L = assemLoad_LFE(Mesh,P7O6(),F_HANDLE);
    [U,FreeDofs] = assemDir_LFE(Mesh,-1,GD_HANDLE);
    L = L - A*U;
    if(~isempty(FreeDofs))
      U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    end
    
    % Real and estimated errors
    
    Err(iter) = H1SErr_LFE(Mesh,U,P7O6(),UEX_HANDLE);
    Mesh = add_Edge2Elem(Mesh);
    Eta_K = ErrEst_HIER(U,Mesh,@Res_Lapl,@STIMA_Lapl_ErrEst,P7O6(),F_HANDLE);  
    Eta(iter) = sqrt(sum(Eta_K.^2));
    
    % Generate plots of error distributions
    
    if(~isempty(find(iter == PLOTS)))
        
        % Generate figure
        
        fig = plot_EST(Eta_K,Mesh);
        colorbar;
        title(['{\bf Error distribution for iteration ' int2str(iter) ' (HIER)}']);
        xlabel(['{\bf' sprintf('Eta = %2.2e',Eta(iter)) '}']);
        text(.1,-.4,'# Coordinates','FontWeight','bold');
        text(.6,-.4,':','FontWeight','bold');
        text(.7,-.4,int2str(size(Mesh.Coordinates,1)),'FontWeight','bold');
        text(.1,-.55,'# Elements','FontWeight','bold');
        text(.6,-.55,':','FontWeight','bold');
        text(.7,-.55,int2str(size(Mesh.Elements,1)),'FontWeight','bold');
        text(.1,-.7,'# Edges','FontWeight','bold');    
        text(.6,-.7,':','FontWeight','bold');
        text(.7,-.7,int2str(size(Mesh.Edges,1)),'FontWeight','bold');  
        
        % Print  to .eps file
        
        FileName = ['ErrEst_HIER_' int2str(iter) '.eps'];
        print('-depsc', FileName);
        close(fig);
        system(['gv ' FileName ' &']);
        
    end

    % Mark elements for refinement
    
    Eta_max = max(Eta_K);
    MarkedElem = find(Eta_K > Eta_max*THETA);
    
    % Refine mesh by largest edge bisection
    
    tmp_Mesh = Mesh;
    Mesh = refine_LEB(Mesh,MarkedElem);
        
    % Update mesh data structure
    
    Mesh = add_Edges(Mesh);
    Mesh = rmfield(Mesh,'BdFlags');
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
    Mesh.BdFlags(Loc) = -1; 
    
    % Print out information
    
    fprintf('Eta(HIER)  :  %2.4e,   iter  :  %d\n',Eta(iter),iter);
  
  end
  
  % Update data
 
  load data.mat;
  Data.nDofs(:,3) = nDofs;
  Data.Eta(:,3) = Eta;
  Data.Err(:,3) = Err;
  Data.LegendName{3} = 'HIER';
  save data.mat Data;
  
  % Clear memory
  
  clear all;
  
  %--------------------------------------------------%
  %        Discrete error for RED refinement         %
  %--------------------------------------------------%
  
  % Initialize constants
  
  load workspace.mat;
  
  % Initialize mesh
  
  Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat'); 
  Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);         
  Mesh = add_Edges(Mesh);                                
  Loc = get_BdEdges(Mesh);
  Mesh.BdFlags = zeros(size(Mesh.Edges,1),1); 
  Mesh.BdFlags(Loc) = -1;
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
    
    Err(iter) = H1SErr_LFE(Mesh,U,P7O6(),UEX_HANDLE);
   
    % Generate new mesh
    
    Mesh = refine_REG(Mesh);  
    
    % Print out information
    
    fprintf('Err (RED)  :  %2.4e,   iter  :  %d\n',Err(iter),iter);
  
  end
  
  % Update data
 
  load data.mat;
  Data.nDofs_RED = nDofs;
  Data.Err_RED = Err;
  save data.mat Data;
  
  % Clear memory
  
  clear all;
  
  % Load data
    
  load data.mat;
  
  % Generate figures
    
  fig = figure();
  plot(Data.nDofs(:,1),Data.Eta(:,1),'-x', ...
       Data.nDofs(:,2),Data.Eta(:,2),'-*', ...
       Data.nDofs(:,3),Data.Eta(:,3),'->', ...
       Data.nDofs_RED,Data.Err_RED,'-o');   
  grid('on');
  set(gca,'XScale','log','YScale','log');
  title('{\bf Convergence rates for energy error and estimators}');
  xlabel('{\bf Dofs [log]}');
  ylabel('{\bf Errors [log]}');
  legend(['Err. est. ' Data.LegendName{1} ' ad. ref.'], ...
         ['Err. est. ' Data.LegendName{2} ' ad. ref.'], ...
         ['Err. est. ' Data.LegendName{3} ' ad. ref.'], ...
         'Disc. err. reg. ref.', ...
         'Location','NorthEast');
  add_Slope(gca,'SouthWest',-0.5);
  add_Slope(gca,'East',-0.3);
  print('-depsc','rate_ErrEst.eps');
  close(fig);
  
  !gv rate_ErrEst.eps &
  
  fig = figure();
  plot(Data.nDofs(:,1),Data.Eta(:,1)./Data.Err(:,1),'-x', ...
       Data.nDofs(:,2),Data.Eta(:,2)./Data.Err(:,2),'-*', ...
       Data.nDofs(:,3),Data.Eta(:,3)./Data.Err(:,3),'->');
  grid('on');
  set(gca,'XScale','log');
  title('{\bf Ratio between energy error and estimators}');
  xlabel('{\bf Dofs [log]}');
  ylabel('{\bf Ratio}');
  legend(['Err. est. ' Data.LegendName{1} ' ad. ref.'], ...
         ['Err. est. ' Data.LegendName{2} ' ad. ref.'], ...
         ['Err. est. ' Data.LegendName{3} ' ad. ref.'], ...
         'Location','East');
  print('-depsc','ratio_ErrEst.eps');
  close(fig);
  
  !gv ratio_ErrEst.eps &
  
  % Delete temp data files
  
  delete *.mat;
  
  % Clear memory
  
  clear all;