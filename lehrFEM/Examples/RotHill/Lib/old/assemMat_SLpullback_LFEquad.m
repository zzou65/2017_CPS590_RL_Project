function varargout= assemMat_SLpullback_LFEquad(Mesh, XtauMesh,intersec,varargin)
% 
% assemble the stiffness matrix of pullback for the semi-lagrange method using one
% point numerical integration on every intersection of triangles of the
% original mesh and the deformed mesh.
%
%  Mesh:           struct to fixed (and straigth) mesh
%  XtauMesh:   struct to  transported mesh
%  intersec:      struct that contains information which element if Mesh
%                       is intersected with which element of XtauMesh
%
%   Copyright 2008-2009 Holger Heumann, Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
tol = 0;

nCoordinates = size(Mesh.Coordinates,1);
nEdges =size(Mesh.Edges,1);
nElements=size(Mesh.Elements,1);

% compute the total number of patches
nPatches = 0;
for i = 1:nElements
    nPatches = nPatches + intersec(i).nElems;
end
 
% Preallocate memory
I = zeros(9*nPatches,1);
J = zeros(9*nPatches,1);
A = zeros(9*nPatches,1);

% Check for element flags
if (isfield(Mesh,'ElemFlag')), flags = Mesh.ElemFlag;
else flags = zeros(nElements,1); end

% Assemble element contributions
loc = 1:9;
% loop over all element of Mesh
for i = 1:nElements
     
   % test whether there is an intersection with some element in XtauMesh,
   if (intersec(i).nElems ~= 0 )
       
       % Vertices of element i
       vid = Mesh.Elements(i,:);
       a1 = Mesh.Coordinates(vid(1),:);
       a2 = Mesh.Coordinates(vid(2),:);
       a3 = Mesh.Coordinates(vid(3),:);
       
       % Compute element mapping
       bK = a1;
       BK = [a2-bK; ...
           a3-bK];
       inv_BK = inv(BK);
       TK = transpose(inv_BK);

        % compute all the polygons of intersections of element i in Mesh with element
        % in XtauMesh 
        polygons = patches(Mesh,XtauMesh,i,intersec(i));
        % loop over all such polygons
        for j = 1:intersec(i).nElems
            
            % Preallocate memory for local stiffness matrix
            Aloc=zeros(3,3);

            % element number of that element in XtauMesh that creates polygon with
            % element i of Mesh
            p_Element = polygons(j).Elem;
            % polygon data
            polygon = polygons(j).Polygon;
           
            % test if area is larger than tolerance tol
            if polyarea(polygon(:,1),polygon(:,2))>tol
                % set barycenter of polygon as evaluation point
                evalpoint = sum(polygon,1)/size(polygon,1);
                
                % evaluate basisfunction of element i at evaluation point
                stevalpoint = (evalpoint - a1)*inv_BK;
                N = shap_LFE(stevalpoint);
                
                % Vertices (id and Coordinates) of element p_Element in XtauMesh
                p_vid = XtauMesh.Elements(p_Element,:);
                p_a1 = XtauMesh.Coordinates(p_vid(1),:);
                p_a2 = XtauMesh.Coordinates(p_vid(2),:);
                p_a3 = XtauMesh.Coordinates(p_vid(3),:);

                % Compute element mapping of element p_Element in XtauMesh
                p_bK= p_a1;
                p_BK = [p_a2-p_bK; ...
                    p_a3-p_bK];
                p_inv_BK = inv(p_BK);
                       
                % evaluate basisfunction element p_Element in XtauMesh at evaluation point
                point_hat=(evalpoint-p_bK)*p_inv_BK;
                p_N =shap_LFE(point_hat);

                % Add contributions to stiffness matrix
                Aloc = N'*p_N;

                Aloc = Aloc*polyarea(polygon(:,1),polygon(:,2));

                I(loc) = set_Rows(vid,3);
                J(loc) = set_Cols(p_vid,3);
                A(loc) = Aloc(:);
                loc = loc+9;
            end                
        end
    end
end

% Assign output arguments
if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
else
    varargout{1} = sparse(I(1:loc(1)-1),J(1:loc(1)-1),A(1:loc(1)-1),nCoordinates,nCoordinates);
end
