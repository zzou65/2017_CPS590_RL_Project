% Run script for the convergence rate for W1F finite element solver on the
% unit square domain.

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    clear Mesh;
    
    NREFS = 5;
    
    MU_HANDLE=@(x,varargin)1;
    U_Handle = @(x,varargin)ones(size(x,1),1);
    F_Handle = @(x,varargin)[sin(pi*x(:,2)) sin(pi*x(:,1))];
    GD_Handle = @(x,varargin)[sin(pi*x(:,2)) sin(pi*x(:,1))];
    UEX_Handle = @(x,varargin)[sin(pi*x(:,2)) sin(pi*x(:,1))];
    
    M_W = zeros(NREFS,1);
    Err_L2 = zeros(NREFS,1);
    
    Mesh.Coordinates = [-1 -1;1 -1;1 1; -1 1];
    Mesh.Elements=[1 2 4;2 3 4];
    
    Mesh = add_Edges(Mesh);
    Mesh = add_Edge2Elem(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc)=-1;
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    Mesh = refine_REG(Mesh);
    
    for i=1:NREFS
        
        Mesh = refine_REG(Mesh);
        M_W(i) = get_MeshWidth(Mesh);

        % Assemble MASS matrix and load vector

        M = assemMat_W1F2nd(Mesh,@MASS_W1F2nd,MU_HANDLE, P3O3());
        L = assemLoad_W1F2nd(Mesh,P7O6(),F_Handle);
        cL = assemCochain_1f2nd(Mesh,F_Handle,gauleg(0,1,10),0);

        %         M = assemMat_W1F(Mesh,@MASS_W1F,MU_HANDLE, P3O3());
        %         L = assemLoad_W1F(Mesh,P7O6(),F_Handle);

        % Incorporate Dirichlet boundary data
        nEdges=size(Mesh.Edges(),1);

        % Solve the system

        U = M\L;

        Err_L2(i) = L2Err_W1F2nd(Mesh,cL,P7O6(),UEX_Handle);
        %Err_L2(i) = L2Err_W1F(Mesh,U(1:nEdges),P7O6(),UEX_Handle)
        
        plot_Norm_W1F(cL(1:nEdges),Mesh);
        
    end
    
    figure
    
    plot(M_W,Err_L2,'-^');
    grid on
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    p = polyfit(log(M_W),log(abs(Err_L2)),1);
    add_Slope(gca,'SouthEast',p(1));
    ylabel('\bf\fontsize{14}Discretization error')
    xlabel('\bf\fontsize{14}h')
    title('\bfConvergence rate for discretization error on the square domain')
    legend('\bf\fontsize{10}W1F:L2 Norm')
    
    % Generate .eps files
    
%     print('-depsc', 'rate_W1F_Sqr.eps');
%     !gv rate_W1F_Sqr.eps &
    
    % Clear memory
    
    %clear all;
    
    