function [errs, mw, msh, U] = exp_1(eps, nRef, plotting, output)
% Runs a solution of the problem given in experiment 1.
% 
% Input arguments:
% eps: epsilon
% nRef: Number of refinements of mesh
% plotting: Plot the solution? (true/false)
% output: Give output to terminal? (true/false)
%
% Output arguments:
% errs: Discretization errors in three norms (L1, L2 and Linf)
% mw: The meshwidth
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
    msh = load_Mesh('mesh_exp_1_coords.dat', 'mesh_exp_1_elements.dat');
  %  msh.Coordinates = [0 0;1 1;1 -1];
  %  msh.Elements = [1 2 3];
    msh.ElemFlag = ones(size(msh.Elements,1),1);
    msh = add_Edges(msh);
    loc = get_BdEdges(msh);
    msh.BdFlags = zeros(size(msh.Edges,1),1);
    msh.BdFlags(loc) = -ones(size(loc));
    for i=1:nRef
        msh = refine_REG(msh);
    end
    msh = add_MidPoints(msh, 'barycentric');

    u = @(x,varargin)(x(1)-(exp(-(1-x(1))/eps)-exp(-1/eps))/(1-exp(-1/eps)));
    k = @(x,varargin)eps;
    c = @(x,varargin)[1 0];
    r = [];
    d = @(x,varargin)(u(x,varargin));
    n = [];
    f = @(x,varargin)1;
    conDom = 1;

    [A, U, L, fd] = assemMat_CD_LFV(msh, k, c, r, d, n, f, conDom, output);

    out('Solving system.');
    U(fd) = A(fd,fd)\L(fd);

    if plotting 
        out('Plotting.');
        plot_LFE(U, msh);
        %plot_Mesh(msh);
        % colorbar;
    end

    out('Calculating errors.');
    errs = [L1Err_LFV(msh, U, P7O6(), u);  L2Err_LFV(msh, U, P7O6(), u); LInfErr_LFV(msh, U, u)];
    out(['Discretization errors: ' num2str(errs(1)) ', ' num2str(errs(2)) ', ' num2str(errs(3)) '.']);

    out('Calculating meshwidth.');
    mw = get_MeshWidth(msh);


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
