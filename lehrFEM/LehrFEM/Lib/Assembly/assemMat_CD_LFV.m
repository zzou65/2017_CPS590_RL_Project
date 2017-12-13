function [A, U, L, fd] = assemMat_CD_LFV(msh, kHandle, cHandle, rHandle, dHandle, nHandle, fHandle, conDom, output)

    if nargin < 9 || isempty(output)
        output = 0;
    end

    tol = 1;
    rFnc = @rFnc1;
    m = [];

    if (nargin < 8 || isempty(conDom)) && ~isempty(cHandle)
        out('Checking for convection domination.');
        
        m = 0;
        for e=1:size(msh.Elements,1)
            vertices = msh.Coordinates(msh.Elements(e,:),:);
            midPoints(1,:) = msh.MidPoints(msh.Vert2Edge(msh.Elements(e,1),msh.Elements(e,2)),:);
            midPoints(2,:) = msh.MidPoints(msh.Vert2Edge(msh.Elements(e,2),msh.Elements(e,3)),:);
            midPoints(3,:) = msh.MidPoints(msh.Vert2Edge(msh.Elements(e,3),msh.Elements(e,1)),:);
            center = msh.CenterPoints(e,:);
            
            g = [abs(dot(novec([midPoints(1,:);center]), cHandle((midPoints(1,:)+center)/2)));
                 abs(dot(novec([midPoints(2,:);center]), cHandle((midPoints(2,:)+center)/2)));
                 abs(dot(novec([midPoints(3,:);center]), cHandle((midPoints(3,:)+center)/2)))];
            d = [norm(vertices(1,:)-vertices(2,:));
                 norm(vertices(2,:)-vertices(3,:));
                 norm(vertices(3,:)-vertices(1,:))];
            u = [kHandle((midPoints(1,:)+center)/2);
                 kHandle((midPoints(2,:)+center)/2);
                 kHandle((midPoints(3,:)+center)/2)];
            z = g.*d./u;
            m = max([m; z]);
        end

        conDom = (m >= tol);
    elseif (nargin < 8 || isempty(conDom)) && isempty(cHandle)
        conDom = 0;
    end

    if conDom
        if isempty(m)
            out(['Problem is locally convection dominated.']);
        else
            out(['Problem is locally convection dominated, m = ' num2str(m) '.']);
        end
    else
        if isempty(m)
            out(['Problem is globally diffusion dominated.']);
        else
            out(['Problem is globally diffusion dominated, m = ' num2str(m) '.']);
        end
    end

    out('Assembling stiffness matrices.');

    A = assemMat_LFV(msh, @STIMA_GenLapl_LFV, kHandle);
    if ~isempty(cHandle)
        %A = A + assemMat_LFV(msh, @STIMA_GenGrad2_LFV, cHandle);
        A = A + assemMat_LFV(msh, @STIMA_GenGrad_LFV, cHandle, kHandle, rFnc, conDom);
    end
    if ~isempty(rHandle)
        A = A + assemMat_LFV(msh, @STIMA_ZerOrd_LFV, rHandle);
    end

    out('Assembling load vector.');

    L = assemLoad_LFV(msh, fHandle);
    if ~isempty(nHandle)
        L = assemNeu_LFV(msh, -2, L, gauleg(0, 1, 4, 1e-6), @(x,varargin)(nHandle(x,varargin)*kHandle(x,varargin)));
    end
    [U, fd] = assemDir_LFV(msh, -1, dHandle);
    L = L - A*U;

    % Output function

    function out(str, level)
        if output
            if nargin < 2 || isempty(level)
                level = 1;
            end
            for i = 1:level
                str = ['- ' str];
            end
            disp(str);
        end
    end

    % R-functions

    function r = rFnc0(z)
        r = 1/2;
    end

    function r = rFnc1(z)
        r = (sign(z) + 1)/2;
    end

    function r = rFnc2(z)
        if z==0, r = 1/2;
        else
            t = max(0, 1-2/abs(z));
            if z<0, r = (1-t)/2;
            else r = (1+t)/2; end
        end
    end

    function r = rFnc3(z)
        if z==0
            r = 1/2;
        else
            r = 1-(1-z/(exp(z)-1))/z;
        end
    end

    % Helping functionss

    function out = novec(v)
        out = v(2,:)-v(1,:);
        if norm(out) ~= 0
            out = [-out(2) out(1)]/norm(out);
        end
    end

end
