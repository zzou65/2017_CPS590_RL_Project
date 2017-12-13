% Run script for convergence rate for the apriori method

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Initialize constant
    
    F_HANDLE = @f_LShap;             % Right hand side source term
    GD_HANDLE = @g_D_LShap;          % Dirichlet boundary data
    UEX_HANDLE = @grad_uex_LShap;    % Exact solution
    FHandle2 = @(x,varagin)x;        % Algebraically graded 2
    BStep = 16;                      % Step for Beta
    Beta = linspace(1,4,BStep);
    MIN_NODES = 5;
    MAX_NODES = 65;
    MIN_NREFS = 3;
    MAX_NREFS = 8;
    
    % Preallocate memory
    
    nDofs = zeros((MAX_NODES-MIN_NODES)/5,BStep);
    Err_AP = zeros((MAX_NODES-MIN_NODES)/5,BStep);
    nDofs_REG = zeros(MAX_NREFS-MIN_NREFS,1);
    Err_REG = zeros(MAX_NREFS-MIN_NREFS,1);
    p = zeros(BStep,2);
    
    % The graded mesh
    
    for iter = 1:BStep
            
        FHandle = @(x,varagin)x.^Beta(iter);   % Algebraically graded
    
        for nNodes = MIN_NODES:5:MAX_NODES
        
            % Generate graded mesh for the reference element
        
            Mesh_REF = graded_RefElem(nNodes,FHandle);
    
            % Generate part1
    
            Mesh_1 = affine_map(Mesh_REF,[ 0 0 ; 1 0 ; 0 1]); 
            Mesh_1.Elements = delaunayn(Mesh_1.Coordinates);
            Mesh_1 = orient_Elems(Mesh_1);
            Mesh_1 = add_Edges(Mesh_1);
    
            Mesh_2 = affine_map(Mesh_REF,[ 0 0 ; 0 1 ; -1 0]);
            Mesh_2.Elements = delaunayn(Mesh_2.Coordinates);
            Mesh_2 = orient_Elems(Mesh_2);
            Mesh_2 = add_Edges(Mesh_2);
  
            Mesh_3 = affine_map(Mesh_REF,[ 0 0 ; -1 0 ; 0 -1]);
            Mesh_3.Elements = delaunayn(Mesh_3.Coordinates);
            Mesh_3 = orient_Elems(Mesh_3);  
            Mesh_3 = add_Edges(Mesh_3);
    
            % Generate part2
    
            Mesh_REF2 = graded_RefElem(nNodes,FHandle2);
    
            Mesh_4 = affine_map(Mesh_REF2,[ 1 1 ; 0 1 ; 1 0]);
            Mesh_4.Elements = delaunayn(Mesh_4.Coordinates);
            Mesh_4 = orient_Elems(Mesh_4);
            Mesh_4 = add_Edges(Mesh_4);
    
            Mesh_5 = affine_map(Mesh_REF2,[ -1 1 ; -1 0 ; 0 1]);
            Mesh_5.Elements = delaunayn(Mesh_5.Coordinates);
            Mesh_5 = orient_Elems(Mesh_5);
            Mesh_5 = add_Edges(Mesh_5);
    
            Mesh_6 = affine_map(Mesh_REF2,[ -1 -1 ; 0 -1; -1 0]);
            Mesh_6.Elements = delaunayn(Mesh_6.Coordinates);
            Mesh_6 = orient_Elems(Mesh_6);
            Mesh_6 = add_Edges(Mesh_6);
    
            % Merge the meshes
    
            Mesh = merge_Mesh(Mesh_1,Mesh_2);
            Mesh = merge_Mesh(Mesh,Mesh_3);    
            Mesh = merge_Mesh(Mesh,Mesh_4);
            Mesh = merge_Mesh(Mesh,Mesh_5);
            Mesh = merge_Mesh(Mesh,Mesh_6); 

            % Error computation
       
            nDofs((nNodes-MIN_NODES)/5+1,iter) = size(Mesh.Coordinates,1);  
      
            % Compute FE solution
  
            Loc = get_BdEdges(Mesh);
            Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
            Mesh.BdFlags(Loc) = -1;
            Mesh.ElemFlag = ones(size(Mesh.Elements)); 

            A = assemMat_LFE(Mesh,@STIMA_Lapl_LFE);
            L = assemLoad_LFE(Mesh,P7O6(),F_HANDLE);
            [U,FreeDofs] = assemDir_LFE(Mesh,-1,GD_HANDLE);
            L = L - A*U;
            if(~isempty(FreeDofs))
                U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
            end
       
            % Real error
    
            Err_AP((nNodes-MIN_NODES)/5+1,iter) = H1SErr_LFE(Mesh,U,P7O6(),UEX_HANDLE);
    
            nNodes
            
        end
        
        clear Mesh Mesh* A L U;
        iter
    end
    
    % The REG refine mesh
    
    % Initialize mesh
  
    Mesh_REG = load_Mesh('Coord_LShap.dat','Elem_LShap.dat'); 
    Mesh_REG.ElemFlag = ones(size(Mesh_REG.Elements,1),1);         
    Mesh_REG = add_Edges(Mesh_REG);                                
    Loc = get_BdEdges(Mesh_REG);
    Mesh_REG.BdFlags = zeros(size(Mesh_REG.Edges,1),1); 
    Mesh_REG.BdFlags(Loc) = -1;
    
    for iter = 1:MAX_NREFS
  
        if(iter>=MIN_NREFS)
        
          nDofs_REG(iter-MIN_NREFS+1) = size(Mesh_REG.Coordinates,1); 
   
          % Compute FE solution
  
          A = assemMat_LFE(Mesh_REG,@STIMA_Lapl_LFE);
          L = assemLoad_LFE(Mesh_REG,P7O6(),F_HANDLE);
          [U,FreeDofs] = assemDir_LFE(Mesh_REG,-1,GD_HANDLE);
          L = L - A*U;
          if(~isempty(FreeDofs))
              U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
          end
    
          % Compute RED discrete errors
    
          Err_REG(iter-MIN_NREFS+1) = H1SErr_LFE(Mesh_REG,U,P7O6(),UEX_HANDLE);
        
        end
   
        % Generate new mesh
    
        Mesh_REG = refine_REG(Mesh_REG);  

    end
        
    % Generate figure
    
    fig = figure();
    plot(nDofs(:,1),Err_AP(:,1),'-x',...
         nDofs(:,2),Err_AP(:,2),'-p',...
         nDofs(:,5),Err_AP(:,5),'-+',...
         nDofs(:,9),Err_AP(:,9),'-*',...
         nDofs(:,12),Err_AP(:,12),'-o',...
         nDofs(:,16),Err_AP(:,16),'-d',...
         nDofs_REG,Err_REG,'-^');   
    grid('on');
    set(gca,'XScale','log','YScale','log');
    title(['{\bf Convergence rate for energy error}']);
    xlabel('{\bf Dofs [log]}');
    ylabel('{\bf Errors [log]}');
    legend(['{\bf\beta = ' num2str(Beta(1),2) '}'],['{\bf\beta = ' num2str(Beta(2),2) '}'],...
           ['{\bf\beta = ' num2str(Beta(5),2) '}'],['{\bf\beta = ' num2str(Beta(9),2) '}'],...
           ['{\bf\beta = ' num2str(Beta(12),2) '}'],['{\bf\beta = ' num2str(Beta(16),2) '}'],...
           'REG', 'Location','NorthEast');
    add_Slope(gca,'SouthWest',-0.5);
    add_Slope(gca,'East',-0.3);

    print('-depsc','rate_GradedErr.eps');
    close(fig);
    
    !gv rate_GradedErr.eps &
    
    for i = 1:BStep
        
        p(i,:) = polyfit(log(nDofs(:,i)),log(Err_AP(:,i)),1)
    
    end
    
    fig = figure();
    plot(Beta,p(:,1),'-x')
    grid('on');
    title('{\bf Relationship between beta and convergence rate}');
    xlabel('{\bf \beta}');
    ylabel('{\bf Convergence rate}');
    print('-depsc','rate_Beta.eps');
    close(fig);
    
    !gv rate_Beta.eps &
       
    % Clear Memory
    
    clear all;