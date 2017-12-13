% Run script for the convergence rate for W1F finite element solver
% including a convection term on the unit square [-1,1]^2

%   Copyright 2009-2009 Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    clear Mesh;
    
    % Initialize constant
    
    NREFS = 6;
    MU_HANDLE=@(x,varargin)1;
    U_Handle = @(x,varargin)ones(size(x,1),1);
    F_Handle = @F;
    W_Handle = @(x,varargin)[2*x(:,1) x(:,2)];
    GD_Handle = @(x,varargin)[sin(pi*x(:,2)) sin(pi*x(:,1))];
    UEX_Handle = @(x,varargin)[sin(pi*x(:,2)) sin(pi*x(:,1))];
    UEX_CURL_Handle = @(x,varargin)-pi*(cos(pi*x(:,1))-cos(pi*x(:,2)));
    
    % Preallocate memory
    
    M_W = zeros(NREFS,1);
    Err_L2 = zeros(NREFS,1);
    Err_HCS = zeros(NREFS,1);
    
    Mesh.Coordinates = [-1 -1;1 -1;1 1;-1 1];
    Mesh.Elements = [1 2 4;2 3 4];
    Mesh = add_Edges(Mesh);
    Mesh = add_Edge2Elem(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;
    Mesh.ElemFlag = ones(size(Mesh.Elements,1),1);
    
    for i=1:NREFS
        
        Mesh = refine_REG(Mesh);
        M_W(i) = get_MeshWidth(Mesh);
        
        [IC,JC,C] = assemMat_W1F(Mesh,@STIMA_Curl_W1F,U_Handle,P7O6());
        [IM,JM,M] = assemMat_W1F(Mesh,@MASS_W1F,MU_HANDLE, P3O3());
        [ICo,JCo,Co] = assemMat_W1F(Mesh,@CONV_Curl_W1F,W_Handle, P3O3());
        A = sparse([IC;IM;ICo],[JC;JM;JCo],[C;M;Co]);
        L = assemLoad_W1F(Mesh,P7O6(),F_Handle);

        % Incorporate Dirichlet boundary data

        [U,FreeDofs] = assemDir_W1F(Mesh,-1,GD_Handle,gauleg(0,1,1));
        L = L - A*U;

        % Solve the system

        U(FreeDofs) = A(FreeDofs,FreeDofs)\L(FreeDofs);
        
        % compute the error
        Err_L2(i) = L2Err_W1F(Mesh,U,P7O6(),UEX_Handle);
        Err_HCS(i) = HCurlSErr_W1F(Mesh,U,P7O6(),UEX_CURL_Handle);
        
    end
    
    figure
    
    plot(M_W,Err_L2,'-^',M_W,Err_HCS,'-*');
    grid on
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    p = polyfit(log(M_W),log(abs(Err_L2)),1);
    add_Slope(gca,'SouthEast',p(1));
    ylabel('\bf\fontsize{14}Discretization error')
    xlabel('\bf\fontsize{14}h')
    title('\bfConvergence rate for discretization error on the disk domain')
    legend('\bf\fontsize{10}W1F:L2 Norm','\bf\fontsize{10}W1F:H(Curl) Semi-Norm')
    
    % Generate .eps files
    
%     print('-depsc', 'rate_W1F_Disk.eps');
%     !gv rate_W1F_Disk.eps &
%     
    % Clear memory
    
    clear all;
    
    
    