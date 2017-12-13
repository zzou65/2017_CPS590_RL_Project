function varargout = assemMat_Contr1f_FEM(Mesh, vHandle, QuadRule,varargin)
% assemMat_Contr1f_FEM assemble contraction of one-forms.
%
%   A = assemMat_Contr1f(Mesh, vHandle, varargin)
%   A = ASSEMMat_Contr1f(MESH,  vHandle) .... and
%   returns the matrix in a sparse representation.
%
%   [I,J,A] = ASSEMMat_Contr1f(MESH, vHandle) .... assembles the global matrix
%   and returns the matrix in an array representation.
%
%
%   Example:
%
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%   V_Handle=@(x,varargin)[ones(size(x,1),1) 0.5.*ones(size(x,1),1)]
%   A = assemMat_Contr1f(Mesh,V_Handle);
%
%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constants
TOL=eps;

nCoordinates = size(Mesh.Coordinates,1);
nEdges =size(Mesh.Edges,1);
nElements=size(Mesh.Elements,1);

% Preallocate memory

I = zeros(9*nElements,1);
J = zeros(9*nElements,1);
A = zeros(9*nElements,1);

% look for upwind nodes on each triangle
loc=1:9;
for i = 1:nElements

    B=zeros(3,3);

    % Vertices
    vid = Mesh.Elements(i,:);
    a1 = Mesh.Coordinates(vid(1),:);
    a2 = Mesh.Coordinates(vid(2),:);
    a3 = Mesh.Coordinates(vid(3),:);

    % Extract global edge numbers

    eidx = [Mesh.Vert2Edge(Mesh.Elements(i,2),Mesh.Elements(i,3)) ...
        Mesh.Vert2Edge(Mesh.Elements(i,3),Mesh.Elements(i,1)) ...
        Mesh.Vert2Edge(Mesh.Elements(i,1),Mesh.Elements(i,2))];

    % Determine the orientation

    if (Mesh.Edges(eidx(1),1)==vid(2)),  p1 = 1;  else   p1 = -1;  end
    if (Mesh.Edges(eidx(2),1)==vid(3)),  p2 = 1;  else   p2 = -1;  end
    if (Mesh.Edges(eidx(3),1)==vid(1)),  p3 = 1;  else   p3 = -1;  end

    % Compute element mapping

    bK = a1;
    BK = [a2-bK; ...
        a3-bK];
    det_BK = abs(det(BK));
    inv_BK = inv(BK);

    % basis functions 0-forms
    lambda = shap_LFE(QuadRule.x);
    
    % basis function 1-forms
    % Compute constant gradients of barycentric coordinate functions
    g1 = [a2(2)-a3(2);a3(1)-a2(1)]/det_BK;
    g2 = [a3(2)-a1(2);a1(1)-a3(1)]/det_BK;
    g3 = [a1(2)-a2(2);a2(1)-a1(1)]/det_BK;

    % Get barycentric coordinates of quadrature points
    nPoints = size(QuadRule.w,1);
    baryc= [QuadRule.x,1-sum(QuadRule.x,2)];

    % Quadrature points in actual element
    % stored as rows of a matrix
    x = QuadRule.x*BK + ones(nPoints,1)*a1;

    % Evaluate coefficient function at quadrature nodes
    Fval = vHandle(x,varargin{:});

    % Evaluate basis functions at quadrature points
    % the rows of b(i) store the value of the the i-th
    % basis function at the quadrature points
    b1 = baryc(:,2)*g3'-baryc(:,3)*g2';
    b2 = baryc(:,3)*g1'-baryc(:,1)*g3';
    b3 = baryc(:,1)*g2'-baryc(:,2)*g1';
    
    w=QuadRule.w*det_BK;
    B(1,1)=p1*sum(w.*lambda(:,1).*sum(b1.*Fval,2));
    B(1,2)=p2*sum(w.*lambda(:,1).*sum(b2.*Fval,2));
    B(1,3)=p3*sum(w.*lambda(:,1).*sum(b3.*Fval,2));
    B(2,1)=p1*sum(w.*lambda(:,2).*sum(b1.*Fval,2));
    B(2,2)=p2*sum(w.*lambda(:,2).*sum(b2.*Fval,2));
    B(2,3)=p3*sum(w.*lambda(:,2).*sum(b3.*Fval,2));
    B(3,1)=p1*sum(w.*lambda(:,3).*sum(b1.*Fval,2));
    B(3,2)=p2*sum(w.*lambda(:,3).*sum(b2.*Fval,2));
    B(3,3)=p3*sum(w.*lambda(:,3).*sum(b3.*Fval,2));

    I(loc) = set_Rows(vid,3);
    J(loc) = set_Cols(eidx,3);
    A(loc) = B(:);
    loc = loc+9;
end


% Assign output arguments

if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
else
    varargout{1} = sparse(I,J,A,nCoordinates,nEdges);
end

return
