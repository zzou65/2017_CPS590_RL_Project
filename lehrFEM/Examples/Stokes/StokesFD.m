function [ u1, u2, p ] = StokesFD( N, f1, f2 )
% ===================================================

% StokesFD
% ========

% Stokes Problem by means of Finite Difference Method

% INPUTS
% ======

% N: number of points in each direction.
% f: right-hand side function handle.

% OUTPUTS
% =======

% u1, u2: the velocity fields components from u = (u1, u2)
% p: the preassure

% ===================================================

% step:
h = 1/N;
x = 0:h:1;

% number of points in each direction:
unk = N - 2;

% number of unknowns: u = (u1, u2) and p;
n = 3 * unk^2;

% "-Lapl" matrix for u1 and u2.
A = gallery('poisson', unk);

% dp/dx1 difference matrix:
e = ones(unk^2, 1);
P1 = spdiags([-h/2*e h/2*e], [-1 1], unk^2, unk^2);

% dp/dx2 difference matrix
P2 = spdiags([-h/2*e h/2*e], [-unk unk], unk^2, unk^2);

% the matrices:
z = sparse(unk^2, unk^2);
H = [A z P1; z A P2; P1' P2' z];

% the right-rand side:
F = zeros(n,1);

for i = 1:size(x)
    for j = 1:size(x)
        F( (j-1)*n+i ) = h^2 * f1( x(i), x(j) );    % first component
        F( (j-1)*n+i + unk^2 ) = h^2 * f2( x(i), x(j) ); % second component
        % third component is zero!
    end
end

% solving the system:
X = H \ F;

u1 = reshape(X(1:unk^2), unk, unk);
u2 = reshape(X(unk^2+1:2*unk^2), unk, unk);
p = reshape(X(2*unk^2+1:end), unk, unk);

% ploting the solution.
fig = figure();
quiver(u1, u2);
print(fig, '-djpeg', sprintf('Velocity field with %d elements', N) );

fig2 = figure();
imagesc(0:unk, 0:unk, p);
print(fig2, '-djpeg', sprintf('Preassure with %d elements', N) );

end