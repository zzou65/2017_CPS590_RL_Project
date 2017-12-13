function L= assemLoad_QFE_SL(Mesh, XtauMesh,intersec,quadflag,QuadRule,FHandle,varargin)
% 
% assemble the load vector (quadratic Lagrangian elements) of pullback for the semi-lagrange method using one
% point numerical integration on every intersection of triangles of the
% original mesh and the deformed mesh.
%
%  Mesh:           struct to fixed (and straigth) mesh
%  XtauMesh:   struct to  transported mesh
%  intersec:      struct that contains information which element if Mesh
%                       is intersected with which element of XtauMesh
%  quadflag:     'bary': use barycenters of polygon 
%                        else : quadrature QuadRule on Delauny triangulaton of mesh 
%  QuadRule:
%
%   Copyright 2008-2011 Holger Heumann, Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
% ASSEMLOAD_QFE Assemble linear FE contributions.
%
%   L = ASSEMLOAD_QFE(MESH,QUADRULE,FHANDLE) assembles the global load 
%   vector for the load data given by the function handle FHANDLE.
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    ELEMFLAG     N-by-1 matrix specifying additional element information.
%
%   QUADRULE is a struct, which specifies the Gauss qaudrature that is used
%   to do the integration:
%    W Weights of the Gauss quadrature.
%    X Abscissae of the Gauss quadrature.
%
%   L = ASSEMLOAD_QFE(MESH,QUADRULE,FHANDLE,FPARAM) also handles the 
%   additional variable length argument list FPARAM to the function handle
%   FHANDLE.
%
%   Example:
%
%   FHandle = @(x,varargin)x(:,1).^2+x(:,2).^2;
%   L = assemLoad_QFE(Mesh,P7O6(),FHandle);
%
%   See also shap_QFE.

%   Copyright 2005-2005 Patrick Meury
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

L = zeros(nCoordinates+nEdges,1);

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
       eidx(1) = Mesh.Vert2Edge(vid(1),vid(2)) + nCoordinates;
       eidx(2) = Mesh.Vert2Edge(vid(2),vid(3)) + nCoordinates;
       eidx(3) = Mesh.Vert2Edge(vid(3),vid(1)) + nCoordinates;
       
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
                       
            % element number of that element in XtauMesh that creates polygon j with
            % element i of Mesh
            p_Element = polygons(j).Elem;
            % boundary vertices of 
            polygon = polygons(j).Polygon;
           
            % test if area is larger than tolerance tol
            if polyarea(polygon(:,1),polygon(:,2))>tol
                
                % test whether we take inexact integration
                if strcmp(quadflag,'bary')

                    % Preallocate memory for local stiffness matrix
                    Aloc=zeros(3,3);

                    % set barycenter of polygon as evaluation point
                    evalpoint = sum(polygon,1)/size(polygon,1);

                    % evaluate basisfunction of element i at evaluation point
                    stevalpoint = (evalpoint - a1)*inv_BK;
                    N = shap_QFE(stevalpoint);
                    
                    % evaluate fhandle at barycenter of element i
                    Fval = FHandle(evalpoint);
                    
                    % Add contributions to stiffness matrix
                    Aloc = N'*Fval*polyarea(polygon(:,1),polygon(:,2));

                    % store entries
                    L(vid) = L(vid)+Aloc(1:3);
                    L(eidx) = L(eidx)+Aloc(4:6);
                    
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
                        
                        % evaluate fhandle at barycenter of element i
                        Fval = FHandle(x);
                        
                        % Add contributions to stiffness matrix
                        for iloc = 1:3
                            L(vid(iloc)) = L(vid(iloc)) + sum(w'.*N(:,iloc).*Fval);
                            L(eidx(iloc)) = L(eidx(iloc)) + sum(w'.*N(:,3+iloc).*Fval);
                        end

                    end  % loop over triangles of polygon
                    
                end % choose quadrature on polygons
            end % if polygon area larger than 0
        end % loop over polygons of element i 
   end % if polygons are not empty 
end % loop over elements of Mesh

return
