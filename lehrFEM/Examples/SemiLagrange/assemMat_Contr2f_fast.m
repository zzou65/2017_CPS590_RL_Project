function varargout = assemMat_Contr2f_new(Mesh, V_Handle, varargin)
% assemMat_Contr2f assemble contraction of two-forms.
% midpointrule for line integration
%
%   A = ASSEMMat_Contr2f(MESH,  v) .... and
%   returns the matrix in a sparse representation.
%
%   [I,J,A] = ASSEMMat_Contr2f(MESH, v) .... assembles the global matrix
%   and returns the matrix in an array representation.
%
%
%   Example:
%
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%   V_Handle=@(x,varargin)[ones(size(x,1),1) 0.5.*ones(size(x,1),1)]
%   v=V_Handle(New_Mesh.Coordinates);
%   A = assemMat_Contr2f(Mesh,v);
%
%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constants
TOL=0;
nCoordinates = size(Mesh.Coordinates,1);
nEdges =size(Mesh.Edges,1);
nElements=size(Mesh.Elements,1);

New_Mesh=Mesh;

% Preallocate memory

I = zeros(3*nElements,1);
J = zeros(3*nElements,1);
A = zeros(3*nElements,1);

% look for upwind nodes on each triangle
loc=1:3;

% rotation to the left
rot=[0 -1; 1 0];

for i = 1:nElements

    B=zeros(1,3);

    % Vertices
    vid = Mesh.Elements(i,:);
    a1 = Mesh.Coordinates(vid(1),:);
    a2 = Mesh.Coordinates(vid(2),:);
    a3 = Mesh.Coordinates(vid(3),:);

    % rotated egdes to left = outward pointing normal

    n1=(a3-a2)*rot;
    n2=(a1-a3)*rot;
    n3=(a2-a1)*rot;

    % Extract global edge numbers

    eidx = [Mesh.Vert2Edge(Mesh.Elements(i,2),Mesh.Elements(i,3)) ...
        Mesh.Vert2Edge(Mesh.Elements(i,3),Mesh.Elements(i,1)) ...
        Mesh.Vert2Edge(Mesh.Elements(i,1),Mesh.Elements(i,2))];

    if (Mesh.Edges(eidx(1),1)==vid(2)),  p1 = 1;  else   p1 = -1;  end
    if (Mesh.Edges(eidx(2),1)==vid(3)),  p2 = 1;  else   p2 = -1;  end
    if (Mesh.Edges(eidx(3),1)==vid(1)),  p3 = 1;  else   p3 = -1;  end

    % Compute element mapping

    bK = a1;
    BK = [a2-bK; ...
        a3-bK];
    det_BK = abs(det(BK));
    inv_BK = inv(BK);

    % Extract nodal vectors
    %
    v1 = V_Handle((a2+a3)/2);
    v2 = V_Handle((a3+a1)/2);
    v3 = V_Handle((a1+a2)/2);

    % Compute decision variables Theta

    B=[v1*n1' v2*n2' v3*n3'];
    B=max(B,TOL).*[p1 p2 p3];
    %     % for first edge
    %     if ((v1)*n1'>TOL)
    %         B(1,1)=p1*v1*n1';
    %     end
    %     % for second edge
    %     if ((v2)*n2'>TOL)
    %         B(1,2)=p2*v2*n2';
    %     end
    %     % for third edge
    %     if ((v3)*n3'>TOL)
    %         B(1,3)=p3*v3*n3';
    %     end
    I(loc) = eidx;
    J(loc) = [i,i,i];
    A(loc) =2*B(:)./(det_BK);
    loc = loc+3;

end

% Assign output arguments

if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
else
    varargout{1} = sparse(I,J,A,nEdges,nElements);
end

return
