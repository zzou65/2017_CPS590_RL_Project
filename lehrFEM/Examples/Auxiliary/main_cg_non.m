% Run script with no preconditioner to slove edge element.

%   Copyright 2005-2006 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
    
    clear Mesh;
    
    % Initialize constant

    NREFS = 5;
    Tal = .5;
    TOL = 1e-12;                                   % Stopping criterion
    MAXIT = 2000;                                  % Maximum number of iterations
    NTime = 5;
    U_Handle = @(x,varargin)ones(size(x,1),1);
    F_Handle = @(x,varargin)(pi^2+Tal)*[sin(pi*x(:,2)) sin(pi*x(:,1))];
    GD_Handle = @(x,varargin)[sin(pi*x(:,2)) sin(pi*x(:,1))];
 
    % Initialize mesh
    
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
    end
    
    % Assemble Curl-curl matrix, MASS matrix and load vector
    
    [ICe,JCe,Ce] = assemMat_W1F(Mesh,@STIMA_Curl_W1F,U_Handle,P7O6());
    [IMe,JMe,Me] = assemMat_W1F(Mesh,@MASS_W1F);
    Se = sparse([ICe;IMe],[JCe;JMe],[Ce;Tal*Me]);
    Le = assemLoad_W1F(Mesh,P7O6(),F_Handle);
    
    [ICn,JCn,Cn] = assemMat_LFE2(Mesh,@STIMA_Curl_LFE2,U_Handle,P7O6());
    [IMn,JMn,Mn] = assemMat_LFE2(Mesh,@MASS_LFE2);
    Sn = sparse([ICn;IMn],[JCn;JMn],[Cn;Tal*Mn]);
    
    % Prolongation matrix
    
    P = convert_LFE2_W1F(Mesh);
    
    % Incorporate Dirichlet boundary data
    
    [Ue,FreeDofs_e] = assemDir_W1F(Mesh,-1,GD_Handle,gauleg(0,1,1));
    Le = Le - Se*Ue;
    [Un,FreeDofs_n] = assemDir_LFE2(Mesh,-1,GD_Handle);
    
    Se = Se(FreeDofs_e,FreeDofs_e);
    Le = Le(FreeDofs_e);
    Sn = Sn(FreeDofs_n,FreeDofs_n);
    P = P(FreeDofs_e,FreeDofs_n);

    % Run cg solver
    
    t = cputime;
    U0 = Le;
    [Ue(FreeDofs_e),flag,relres,iter,resvec] = cg_solve(U0,Se,Le,TOL,MAXIT);
    fprintf('Runtime of CG solver [s]  :  %f\n',cputime-t);

    % Plot figures
    
    fig = figure('Name','Non-Preconditioned CG solver');
    plot(resvec,'rx');
    title('{\bfNon-Preconditioned CG solver}');
    xlabel('{\bf Iteration number}');
    ylabel('{\bf Relative residual [log]}');
    set(gca,'YScale','log');
    
    % Output .eps files
    
    print('-depsc', 'CG_non_prec.eps');
    !gv CG_non_prec.eps &
    