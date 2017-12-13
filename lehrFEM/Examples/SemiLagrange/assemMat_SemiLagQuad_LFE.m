function varargout= assemMat_SemiLagQuad_LFE(Mesh,tracedvertices, varargin)
% assemMat_SemiLagQuad_LFE calculates the stiffness matrix for 
% Semi-Lagrange for linear finite elements  using simple vertex based
% quadrature rule
%
% Mesh structure; tracedvertices array with coodinates of traced vertices and
% local local elements.
%
% Example:
% M_0 =  assemMat_Mass0fD(NewMesh);
% pulled_back_vertices=trace_vertices(NewMesh,-dt*directions);
% A =assemMat_SemiLag_LFE(Mesh,pulled_back_vertices);
% B =assemMat_SemiLag_Quad_LFE(Mesh,pulled_back_vertices); 
% % M_0 * A == B;
%
% Copyright 2008-2009 Holger Heumann
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

nCoordinates = size(Mesh.Coordinates,1);
nEdges =size(Mesh.Edges,1);
nElements=size(Mesh.Elements,1);

% Preallocate memory
nElements = size(Mesh.Elements,1);

% Preallocate memory

I = zeros(9*nCoordinates,1);
J = zeros(9*nCoordinates,1);
A = zeros(9*nCoordinates,1);

% Check for element flags
if (isfield(Mesh,'ElemFlag')), flags = Mesh.ElemFlag;
else flags = zeros(nElements,1); end

% Assemble element contributions
loc = 1:9;
for i = 1:nElements
    Aloc=zeros(3,3);

    % Vertices
    vid = Mesh.Elements(i,:);
    a1 = Mesh.Coordinates(vid(1),:);
    a2 = Mesh.Coordinates(vid(2),:);
    a3 = Mesh.Coordinates(vid(3),:);

    % Compute element mapping
    bK = a1;
    BK = [a2-bK; ...
        a3-bK];
    det_BK = abs(det(BK));
    inv_BK = inv(BK);

    % loop over all three vertices
    Map_vid_vidx=zeros(3,3);
    Map_vid_values=zeros(3,3);
    
    N=shap_LFE([0 0; 1 0; 0 1]);
    
    for j=1:3
        point=tracedvertices(vid(j),[1,2]);
        p_Element=tracedvertices(vid(j),3);
        
        % Vertices
        p_vid = Mesh.Elements(p_Element,:);
        p_a1 = Mesh.Coordinates(p_vid(1),:);
        p_a2 = Mesh.Coordinates(p_vid(2),:);
        p_a3 = Mesh.Coordinates(p_vid(3),:);
        
        % Compute element mapping
        p_bK= p_a1;
        p_BK = [p_a2-p_bK; ...
                    p_a3-p_bK];
        p_det_BK = abs(det(p_BK));
        p_inv_BK = inv(p_BK);

        point_hat=(point-p_bK)*p_inv_BK;

        Map_vid_vidx(j,:)=p_vid;
        Map_vid_values(j,:)=shap_LFE(point_hat);
    end
    
    % Add contributions to stiffness matrix

    Aloc = det_BK/6*Map_vid_values; 
    %Aloc = Map_vid_values; 
    I(loc) = set_Rows(vid,3);
    J(loc) = Map_vid_vidx(:);
    A(loc) = Aloc(:);
    loc = loc+9;
end

% Assign output arguments

if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
else
    varargout{1} = sparse(I,J,A,nCoordinates,nCoordinates);
end