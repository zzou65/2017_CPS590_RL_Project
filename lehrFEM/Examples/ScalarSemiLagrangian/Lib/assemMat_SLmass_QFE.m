function varargout= assemMat_SLmass_QFE(Mesh, XtauMesh,intersec,quadflag,QuadRule,varargin)
% 
% assemble the mass matrix (quadratic Lagrangian elements) for the semi-lagrange method using either one
% point numerical integration on every intersection of triangles of the
% original mesh and the deformed mesh or quadrature on some triangulation
% of the resulting polygons. If we use a quadrature that is inexact quadrature 
%
%  Mesh:           struct to fixed (and straigth) mesh
%  XtauMesh:   struct to  transported mesh
%  intersec:      struct that contains information which element if Mesh
%                       is intersected with which element of XtauMesh
%  quadflag:     'bary': use barycenters of polygon 
%                        else : quadrature QuadRule on Delauny triangulaton of mesh 
%  QuadRule:   quadrature rule for triangulation of polygons.
%
%   Copyright 2008-2011 Holger Heumann, Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
tol = 0;

nCoordinates = size(Mesh.Coordinates,1);
nEdges =size(Mesh.Edges,1);
nElements=size(Mesh.Elements,1);
nPts = size(QuadRule.w,1);

% compute the total number of patches
nPatches = 0;
for i = 1:nElements
    nPatches = nPatches + intersec(i).nElems;
end
 
% Preallocate memory
I = zeros(36*nPatches,1);
J = zeros(36*nPatches,1);
A = zeros(36*nPatches,1);

% Check for element flags
if (isfield(Mesh,'ElemFlag')), flags = Mesh.ElemFlag;
else flags = zeros(nElements,1); end

% Assemble element contributions
loc = 1:36;
% loop over all element of Mesh
for i = 1:nElements
     
   % test whether there is an intersection with some element in XtauMesh,
   if (intersec(i).nElems ~= 0 )
       
       % Vertices of element i
       vid = Mesh.Elements(i,:);
       a1 = Mesh.Coordinates(vid(1),:);
       a2 = Mesh.Coordinates(vid(2),:);
       a3 = Mesh.Coordinates(vid(3),:);
       
       % global ids of vertex and edge dofs of element i
       idx = [vid ...
           Mesh.Vert2Edge(Mesh.Elements(i,1),Mesh.Elements(i,2))+nCoordinates ...
           Mesh.Vert2Edge(Mesh.Elements(i,2),Mesh.Elements(i,3))+nCoordinates ...
           Mesh.Vert2Edge(Mesh.Elements(i,3),Mesh.Elements(i,1))+nCoordinates];
       
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
            Aloc=zeros(6,6);

            % element number of that element in XtauMesh that creates polygon j with
            % element i of Mesh
            p_Element = polygons(j).Elem;
            % boundary vertices of polygon
            polygon = polygons(j).Polygon;
           
            % test if area is larger than tolerance tol
            if polyarea(polygon(:,1),polygon(:,2))>tol

                % test whether we take inexact integration
                if strcmp(quadflag,'bary')
                    % set barycenter of polygon as evaluation point
                    evalpoint = sum(polygon,1)/size(polygon,1);

                    % evaluate basisfunction of element i at evaluation point
                    stevalpoint = (evalpoint - a1)*inv_BK;
                    N = shap_QFE(stevalpoint);

                    % Add contributions to stiffness matrix
                    Aloc = N'*N*polyarea(polygon(:,1),polygon(:,2));

                    % store entries
                    I(loc) = set_Rows(idx,6);
                    J(loc) = set_Cols(idx,6);
                    A(loc) = Aloc(:);
                    loc = loc+36;
                    
                else
                    % we use exact integration
                    % Delauny triangulation TODO (meshquality)
                    TRI = delaunay(polygon(:,1),polygon(:,2));
                    nTRI = size(TRI,1);
                    % loop over all elements
                    for k =1:nTRI

                        % Compute element mapping on local element
                        TbK = polygon(TRI(k,1),:);
                        TBK = [polygon(TRI(k,2),:)-TbK; ...
                            polygon(TRI(k,3),:)-TbK];
                        Tdet_BK = abs(det(TBK));

                        % transform evaluation points and weights of quadrature rule
                        x = QuadRule.x*TBK+ones(nPts,1)*TbK;
                        w = QuadRule.w'*Tdet_BK;

                        %evaluate basisfunction of element i at evaluationpoints
                        stevalpoints = (x - ones(nPts,1)*a1)*inv_BK;
                        N = shap_QFE(stevalpoints);

                        % Add contributions to stiffness matrix
                        for iloc = 1:6
                            for jloc = 1:6
                                Aloc(iloc,jloc) = Aloc(iloc,jloc) + sum(w'.*N(:,iloc).*N(:,jloc));
                            end
                        end

                    end  % loop over triangles of polygon
                    
                    % store entries
                    I(loc) = set_Rows(idx,6);
                    J(loc) = set_Cols(idx,6);
                    A(loc) = Aloc(:);
                    loc = loc+36;
                end % choose quadrature on polygons
            end % if polygon area larger than 0
        end % loop over polygons of element i  
   end % if polygons are not empty 
end % loop over elements of Mesh

% Assign output arguments
if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
else
    varargout{1} = sparse(I(1:loc(1)-1),J(1:loc(1)-1),A(1:loc(1)-1),nCoordinates+nEdges,nCoordinates+nEdges);
end