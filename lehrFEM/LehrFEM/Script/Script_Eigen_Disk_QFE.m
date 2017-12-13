% Run script for the Error of eigenvalue problem in the unit disk domain
% Notice: If you want to handle more(e.g. N) eigenvalues, please read the 
%         following:
%         The analytical eigenvalue is the zero point of Bessel function.
%         See http://mathworld.wolfram.com/BesselFunctionZeros.html to 
%         choose the the first to Nth smallest bessel zero point to fill in 
%         the vector ZBessel.
%         And the vector Select should be the smallest Nth from the 
%         following vector:
%         [1 2 4 6 7 9 11 13 15 16 18 20 22 24 26 28 30 31 33 35 37]

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland 

    % Initialize constant
    
    GD_HANDLE = @(x,varargin)zeros(size(x,1),1);
    Select = [1 2 4 6];
    NEigen = 6;
    ZBessel = [2.40482555769577 3.83170597020751 5.13562230184068 5.52007811028631];
    ZBessel = sort(ZBessel);
    Lambda = ZBessel.^2;
    H0 =[ .25 .2 .1 .05 .02 .01 ]';
    NRef = size(H0,1);
    
    % Preallocate memory
    
    M_W = zeros(NRef,1);
    D_LambdaH = zeros(NRef,4);
    
    % Run the script
    
    for i = 1:NRef 
    
        % Initialize mesh
    
        C = [0 0];             % Center of circle
        R = 1;                 % Radius of circle 
        BBOX = [-1 -1; 1 1];   % Bounding box
        DHANDLE = @dist_circ;  % Signed distance function
        HHANDLE = @h_uniform;  % Element size function
        FIXEDPOS = [];         % Fixed boundary vertices of the mesh
        DISP = 0;              % Display flag

        Mesh = init_Mesh(BBOX,H0(i),DHANDLE,HHANDLE,FIXEDPOS,DISP,C,R);  
        Mesh = add_Edges(Mesh);
        Loc = get_BdEdges(Mesh);
        Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
        Mesh.BdFlags(Loc) = -1;
        Mesh.ElemFlag = zeros(size(Mesh.Elements,1),1);

        % Get mesh width
        
        M_W(i) = get_MeshWidth(Mesh);
   
        % Assemble stiffness matrix and mass matrix

        A = assemMat_QFE(Mesh,@STIMA_Lapl_QFE,P7O6());
        M = assemMat_QFE(Mesh,@MASS_QFE,P7O6());


        % Incorporate Dirichlet boundary data

        [U,FreeNodes] = assemDir_QFE(Mesh,-1,GD_HANDLE);
        A = A(FreeNodes,FreeNodes);
        M = M(FreeNodes,FreeNodes);

        % Solve eigenvalue problem

        [ V d ] = eigs(A,M,NEigen,'sm');

        for j = 1:NEigen

            U(FreeNodes,j) = V(:,j);
            norm(A*V(:,j)-d(j,j)*M*V(:,j))

        end

        % Plot the solution

        if ( i == NRef )

            plot_QFE(U(:,NEigen),Mesh);
            colorbar;
            title(['{\bf Solution: 1 '  '  eigenvalue = '...
                num2str(d(NEigen,NEigen))  ', meshwidth = ' num2str(M_W(i)) '}' ]);
            print('-depsc', 'Eig_func1_disk_QFE.eps');
            !gv Eig_func1_disk_QFE.eps &

            plot_QFE(U(:,NEigen-1),Mesh);
            colorbar;
            title(['{\bf Solution: 2 '  '  eigenvalue = '...
                num2str(d(NEigen-1,NEigen-1))  ', meshwidth = ' num2str(M_W(i)) '}' ]);
            print('-depsc', 'Eig_func2_disk_QFE.eps');
            !gv Eig_func2_disk_QFE.eps &

        end

        % Compute the error
        
        d = diag(d);
        d = d(6:-1:1);
        D_LambdaH(i,:) = d(Select)'-Lambda;
        
    end
    
    % Output the .eps file
    
    figure;
    plot(M_W,abs(D_LambdaH(:,1)),'-x',...
         M_W,abs(D_LambdaH(:,2)),'-^',...
         M_W,abs(D_LambdaH(:,3)),'-h',...
         M_W,abs(D_LambdaH(:,4)),'-o'); 
    grid on
    set(gca,'XScale','log','YScale','log','XDir','reverse');
    title('\bfConvergence rate for eigenvalue errors(Disk) QFE');
    xlabel('\bfh');
    ylabel('\bf|\lambda^{h}_{l}-\lambda_{l}|');
    legend('EigenVal.1','EigenVal.2','EigenVal.3','EigenVal.4',...
            'Location','NorthEast')
    p = polyfit(log(M_W),log(abs(D_LambdaH(:,4))),1);
    add_Slope(gca,'SouthEast',p(1));
    print('-depsc', 'rate_Eig_disk_QFE.eps');
    !gv rate_Eig_disk_QFE.eps &
  
    % Clear memory
    
%     clear all;
    
    