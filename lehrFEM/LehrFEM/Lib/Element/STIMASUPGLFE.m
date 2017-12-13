function Aloc = STIMASUPGLFE(Vertices, flag, QuadRule, VHandle, a,d1,d2, varargin)
%   ALOC = STIMA\_SUPG\_LFE(VERTICES) provides the extra terms for SUPG stabilization to be
%  added to the Galerkin element matrix for linear finite elements
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
%  
%   a: diffusivity
%   d1 d2: apriori chosen constants for SUPG-modification
%
%   Flag not used, needed for interface to assemMat\_LFE
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the numerical integration:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.e: 
%
%   VHANDLE is function handle for velocity field   

% Preallocate memory for element matrix
Aloc = zeros(3,3);

% Analytic computation of entries of element matrix using barycentric
% coordinates, see Sect.~\ref{sec:galmatcomp}
l1x = Vertices(2,2)-Vertices(3,2); 
l1y = Vertices(3,1)-Vertices(2,1);
l2x = Vertices(3,2)-Vertices(1,2);
l2y = Vertices(1,1)-Vertices(3,1); 
l3x = Vertices(1,2)-Vertices(2,2);
l3y = Vertices(2,1)-Vertices(1,1);

% Compute element mapping

P1 = Vertices(1,:);
P2 = Vertices(2,:);
P3 = Vertices(3,:);
  
BK = [ P2 - P1 ; P3 - P1 ];          % transpose of transformation matrix
det_BK = abs(det(BK));               % twice the area of the triagle
    
nPoints = size(QuadRule.w,1);
  
% Quadrature points in actual element stored as rows of a matrix
x = QuadRule.x*BK + ones(nPoints,1)*P1;

% Evaluate coefficient function (velocity) at quadrature nodes
c =VHandle(x,varargin{:});
% Entries of anisotropic diffusion tensor
FHandle=[c(:,1).*c(:,1) c(:,1).*c(:,2) c(:,2).*c(:,1) c(:,2).*c(:,2)];

% Compute local PecletNumber for SUPG control parameter
hK=max([norm(P2-P1),norm(P3-P1),norm(P2-P3)]);
v_infK=max(abs(c(:))); PK=v_infK*hK/(2*a);
% Apply quadrature rule and fix constant part
w = QuadRule.w; e = sum((FHandle.*[w w w w]), 1);
te = (reshape(e,2,2)')/det_BK;

% Compute Aloc values
Aloc(1,1) = (te*[l1x l1y]')'*[l1x l1y]';            
Aloc(1,2) = (te*[l1x l1y]')'*[l2x l2y]';            
Aloc(1,3) = (te*[l1x l1y]')'*[l3x l3y]';            
Aloc(2,2) = (te*[l2x l2y]')'*[l2x l2y]';            
Aloc(2,3) = (te*[l2x l2y]')'*[l3x l3y]';            
Aloc(3,3) = (te*[l3x l3y]')'*[l3x l3y]';            
Aloc(2,1) = (te*[l2x l2y]')'*[l1x l1y]';            
Aloc(3,1) = (te*[l3x l3y]')'*[l1x l1y]';            
Aloc(3,2) = (te*[l3x l3y]')'*[l2x l2y]';

if (PK<=1), Aloc=d1*hK^2/a*Aloc;
else Aloc=d2*hK*Aloc; end

return