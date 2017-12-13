function [errs, mw, msh, U] = exp_3(a, nRef, plotting, output)
% Run a solution of the problem given in experiment 3.
%
% Input arguments:
% a: a
% nRef: Number of refinements of the mesh
% plotting: Plot the solution (true/false)
% output: Give output to the terminal (true/false)
%
% Output arguments:
% errs: Vector with discretization errors in three norms (L1, L2, Linf)
% mw: Meshwidth
% msh: The mesh
% U: Finite volume solution vector

    if nargin < 2 || isempty(nRef)
        nRef = 0;
    end
    if nargin < 3 || isempty(plotting)
        plotting = 0;
    end
    if nargin < 4 || isempty(output)
        output = 0;
    end

    out('Generating mesh.');
    msh = load_Mesh('mesh_exp_2_coords.dat', 'mesh_exp_2_elements.dat');
    msh.ElemFlag = ones(size(msh.Elements,1),1);
    msh = add_Edges(msh);
    loc = get_BdEdges(msh);
    msh.BdFlags = zeros(size(msh.Edges,1),1);
    msh.BdFlags(loc) = -ones(size(loc));
    for i=1:nRef
        msh = refine_REG(msh);
    end
    msh = add_MidPoints(msh, 'barycentric');

    u = @ufunc;
    k = @(x,varargin)1;
    c = @cfunc;
    r = [];
    d = @ufunc;
    n = [];
    f = @(x,varargin)0;
    condom = 1;

    [A, U, L, fd] = assemMat_CD_LFV(msh, k, c, r, d, n, f, condom, output);

    out('Solving system.');
    U(fd) = A(fd,fd)\L(fd);

    if plotting 
        out('Plotting.');
        plot_LFE2(U, msh);
        %hold on
        %U_exc = zeros(size(U));
        %for i=1:length(U_exc)
        %    U_exc(i) = u(msh.Coordinates(i,:));
        %end
        %plot_LFE2(U_exc, msh);
        %plot_Mesh(msh);
        % colorbar;
    end

    out('Calculating errors.');
    errs = [L1Err_LFV(msh, U, P7O6(), u);  L2Err_LFV(msh, U, P7O6(), u); LInfErr_LFV(msh, U, u)];
    out(['Discretization errors: ' num2str(errs(1)) ', ' num2str(errs(2)) ', ' num2str(errs(3)) '.']);

    out('Calculating meshwidth.');
    mw = get_MeshWidth(msh);

    % Psi-related functions, very important...
    
    function out = psi(r)
        out = 1/(1+exp(a*(r-1)));
    end

    function out = cfunc(x, varargin)
        r = norm(x);
        if x(1)==0 && x(2)==0
            theta = pi/4;
        elseif x(1)==0
            theta = pi/2;
        else
            theta = atan(x(2)/x(1));
        end
        p = psi(r);
        pg = a*p*(p-1);
        out = [pg*cos(theta) pg*sin(theta)];
    end

    function out = ufunc(x, varargin)
        out = exp(psi(norm(x)));
    end

    % Helping functions

    function out(text, level)
        if output
            if nargin < 2 || isempty(level),
                level = 1;
            end
            for i=1:level
                text = ['- ' text];
            end
            disp([text]);
        end
    end
end
