function varargout = cg(A,b,x0,TOL,MAXIT)

% CG solves a linear complementary problem(LCP)
%
%   X = CG(A,B,X0,TOL,MAXIT) soves the LCP problem given by
%   
%   x'(Ax-b) = 0, x >= 0, Ax-b >= 0 using a conjugate gradient method
% 
%   The tolerence TOL is a stopping criterion for the number of iterations 
%   MAXIT specifies the maximum number of iterations
%
%   [X CONV] = CG(A,B,X0,TOL,MAXIT) returns a flag that indicates weather the
%   iteration converges or not within the specified MAXIT number of
%   iterations
%
%   [X CONV ERR] = CG(A,B,X0,TOL,MAXIT) returns a vector with the relative 
%   residual ERR at each iteration
%
%   [X CONV ERR IT] = CG(A,B,X0,TOL,MAXIT) also returns the number of iterations
%   that were performed
%
%   [X CONV ERR IT IOUT] = CG(A,B,X0,TOL,MAXIT) also returns a vector with
%   the iteration numbers when the outer loop restarts
%
%   Example:
%
%   x = cg(A,b,xo,1e-5,500)

%   Copyright 2006-2006 Kari Borset
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland  



% set parameters
lb   = 1e-10;
n    = length(x0);
conv = 1;
r0 = norm(b);

% initialize indices
I  = (1:n)';
J  = zeros(n,1);
ni = n;
nj = 0;
it = 0;
ERR = [];
iout = [];

% initialixe x
x = x0;

%--------------------------------------------------------------------------
% Outer Iteration
%--------------------------------------------------------------------------
for k=1:MAXIT,
    iout = [iout it];
    y = A*x - b;

    % loop over points to get set I
    ii=0;
    jj=0;
    I0 = I;
    I  = zeros(n,1);
    J  = zeros(n,1);
    for l=1:n,
        if ((abs(x(l))<=lb) & ((abs(y(l))<=lb | y(l) > 0))),
            ii = ii + 1;
            I(ii) = l;
        else
            jj = jj + 1;
            J(jj) = l;
        end
    end

    % Stop if I don't change
    l = find(I-I0);
    if (length(l) == 0),
        break
    else
        ni = ii;
        nj = jj;
    end

    %--------------------------------------------------------------------------
    % Inner Iteration
    %--------------------------------------------------------------------------
    xi  = x(I(1:ni));
    xj  = x(J(1:nj));
    bj  = b(J(1:nj));
    Ajj = A(J(1:nj),J(1:nj));
    Aji = A(J(1:nj),I(1:ni));

    % compute residual
    r = bj - Aji*xi - Ajj*xj;
    p = r;
    for q=1:MAXIT,
        it = it + 1;
        % maximal step size
        s = Ajj*p;
        t = r'*r;
        l = find(p<0);
        a = min([-xj(l)./p(l) ; t/(p'*s)]);
  
        % new values
        xj = xj + a*p;
        r = r - a*s;
        %err = max(abs(r));
        err = norm(r)/r0;
        
        ERR = [ERR err];
        if(err<=TOL)
            x(J(1:nj)) = xj;
            break
        end
        
        % check termination conditions
        ii = 0;
        jj = 0;
        J0 = J;
        I0 = I;
        for l=1:nj,
            if (abs(xj(l))<=lb),
                ii = ii + 1;
                I(ni + ii) = l;
            else
                jj = jj + 1;
                J(jj) = J0(l);
            end
        end
        if (ii>0),
            x(J0(1:nj)) = xj;
            if (ii==nj),
                J = J0;
                I = I0;
                break
            else
                ni = ii;
                nj = jj;
                xi  = x(I(1:ni));
                xj  = x(J(1:nj));
                bj  = b(J(1:nj));
                Ajj = A(J(1:nj),J(1:nj));
                Aji = A(J(1:nj),I(1:ni));

                % compute residual
                r = bj - Aji*xi - Ajj*xj;
                p = r;
                continue
            end
        end
        p = r + (r'*r/t)*p;
    end
end

if k==MAXIT || q==MAXIT,
    conv = 0;
end

switch (nargout)
    case 1
        varargout{1} = x;
    case 2
        varargout{1} = x;
        varargout{2} = conv;
    case 3 
        varargout{1} = x;
        varargout{2} = conv;
        varargout{3} = ERR;
    case 4
        varargout{1} = x;
        varargout{2} = conv;
        varargout{3} = ERR;
        varargout{4} = it;
    case 5
        varargout{1} = x;
        varargout{2} = conv;
        varargout{3} = ERR;
        varargout{4} = it;
        varargout{5} = iout;
end
return