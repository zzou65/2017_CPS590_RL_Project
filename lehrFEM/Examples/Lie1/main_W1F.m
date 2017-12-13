% Run script for W1F finite element solver.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    clear Mesh;
    
    % Initialize constant
    
    NREFS =5;
    EPSI_Handle = @(x,varargin)ones(size(x,1),1);
    FV_Handle = @(x,varargin)[2*ones(size(x,1),1) -2*ones(size(x,1),1)];
    %FV_Handle = @(x,varargin)pi*[cos(pi*x(:,1))-cos(pi*x(:,2)) -cos(pi*x(:,1))+cos(pi*x(:,2))];
    FC_Handle = @(x,varargin)pi^2*[sin(pi*x(:,2)) sin(pi*x(:,1))];
    FM_Handle = @(x,varargin)[sin(pi*x(:,2)) sin(pi*x(:,1))];
    F_Handle=FC_Handle;
    GD_Handle = @(x,varargin)[sin(pi*x(:,2)) sin(pi*x(:,1))];
  % UEX_Handle = @(x,varargin)[sin(pi*x(:,2)) sin(pi*x(:,1))];
    UEX_Handle = @(x,varargin)[-x(:,2) x(:,1)];
    UEX_CURL_Handle = @(x,varargin)pi*cos(pi*x(:,1))-pi*cos(pi*x(:,2));
    MU_HANDLE=@(x,varargin)1;
    V_Handle=@(x,varargin)ones(size(x,1),2) ;
    ZeroHandle=@(x,varargin)zeros(size(x,1),2);
    % Initialize mesh
    
    Mesh.Coordinates = [-1 -1;1 -1;1 1;-1 1];
    Mesh.Elements = [1 2 4;2 3 4];
    Mesh = add_Edges(Mesh);
    Mesh = add_Edge2Elem(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = [-1 -1 -1 -1];
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    Mesh=add_UPWINDDATA(Mesh, V_Handle);
    
    %plot_Mesh(Mesh,'petas');
    for i=1:NREFS
        
%     plot_Mesh(Mesh,'petaso');
%        Start=Mesh.Coordinates(Mesh.Edges(:,1),:);
%        End=Mesh.Coordinates(Mesh.Edges(:,2),:);
%        hold on;
%        quiver(Start(:,1),Start(:,2),End(:,1)-Start(:,1),End(:,2)-Start(:,2),0);
%        hold off;
%     
    
     % Assemble Curl-curl matrix, convection matrix, MASS matrix and load vector
 
 %    C = assemMat_W1F(Mesh,@STIMA_Curl_W1F,EPSI_Handle,P7O6());
     M = assemMat_W1F(Mesh,@MASS_W1F,MU_HANDLE, P7O6());
     V = assemMat_W1F(Mesh,@Contract_1F,V_Handle);
     
     Mesh=add_UPWINDDATA(Mesh, V_Handle);
     V2 = assemMat_Lie_W1F(Mesh,@Contract_1Ffem,V_Handle);
     
     U_L = assemLoad_W1F(Mesh,P7O6(),UEX_Handle);
     U_e = M\U_L;
     
     F_L = assemLoad_W1F(Mesh,P7O6(),FV_Handle);
     F_e=M\F_L;
  %   plot_W1F(F_e,Mesh);
     [U,FreeDofs] = assemDir_W1F(Mesh,[-1],UEX_Handle,gauleg(0,1,4));
    
     
     e=(U_e'*V2*U_e);
     loc=(get_BdEdges(Mesh));
     nEdges=size(Mesh.Edges(),1);
     %max(e(setdiff(1:nEdges,loc)))
     max(e)
     Mesh = refine_REG(Mesh);
    end
    % Incorporate Dirichlet boundary data
    
    [U,FreeDofs] = assemDir_W1F(Mesh,[-1],GD_Handle,gauleg(0,1,4));
   % plot_W1F(U,Mesh);
    L = L - A*U;
    
    % Solve the system
    
    U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
    fprintf('Runtime of direct solver [s]  :  %f\n',cputime-t);
    % Plot the solution
    
%     plot_Mesh(Mesh,'as')
%      hold on
    plot_W1F(U,Mesh);
    
    % Clear memory
    
 clear all
    
    