function varargout= assemMat_SL_PDG(Mesh, XtauMesh,nDofs,intersec,quadflag,QuadRule, Shap, varargin)
% 
% assemble the pullback and mass matrix for the semi-lagrange method using either one
% point numerical integration on every intersection of triangles of the
% original mesh and the deformed mesh or quadrature on some triangulation
% of the resulting polygons. If we use a quadrature that is inexact quadrature 
%
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
 
% Preallocate memory
% compute the total number of patches
nPatches = 0;
for i = 1:nElements
    nPatches = nPatches + intersec(i).nElems;
end
 
nEntries = nPatches * nDofs^2;

% Preallocate memoryn       
IA = zeros(nEntries,1);
JA = zeros(nEntries,1);
A = zeros(nEntries,1);

IM = zeros(nEntries,1);
JM = zeros(nEntries,1);
M = zeros(nEntries,1);

% Check for element flags
if (isfield(Mesh,'ElemFlag')), flags = Mesh.ElemFlag;
else flags = zeros(nElements,1); end

% Assemble element contributions
%loc = 1:36;
offset = 0;
offsetM = 0 ;
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
       idx = nDofs*(i-1)+(1:nDofs);
       
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
            % boundary vertices of polygon
            polygon = polygons(j).Polygon;
           
            % test if area is larger than tolerance tol
            if polyarea(polygon(:,1),polygon(:,2))>tol
                
                % Vertices (id and Coordinates) of element p_Element in XtauMesh
                p_vid = XtauMesh.Elements(p_Element,:);
                p_a1 = XtauMesh.Coordinates(p_vid(1),:);
                p_a2 = XtauMesh.Coordinates(p_vid(2),:);
                p_a3 = XtauMesh.Coordinates(p_vid(3),:);
                
                % global ids of vertex and edge dofs of element p_Element
                p_idx = nDofs*(p_Element-1)+(1:nDofs);
       
                % Compute element mapping of element p_Element in XtauMesh
                p_bK= p_a1;
                p_BK = [p_a2-p_bK; ...
                    p_a3-p_bK];
                p_inv_BK = inv(p_BK);

                % test whether we take inexact integration
                if strcmp(quadflag,'bary')
                    
                    % set barycenter of polygon as evaluation point
                    evalpoint = sum(polygon,1)/size(polygon,1);

                    % evaluate basisfunction of element i at evaluation point
                    stevalpoint = (evalpoint - a1)*inv_BK;
                    N = Shap(stevalpoint);

                    % evaluate basisfunction element p_Element in XtauMesh at evaluation point
                    point_hat=(evalpoint-p_bK)*p_inv_BK;
                    p_N =shap_hp(point_hat,pmax);
                 
                    % quadrature weight
                    weights =  polyarea(polygon(:,1),polygon(:,2));
                    
                    % Add contributions to stiffness matrix
                    Aloc = N'*p_N*polyarea(polygon(:,1),polygon(:,2));
                    Mloc = N'*N*polyarea(polygon(:,1),polygon(:,2));

                    % store entries
                    loc = offset + (1:nDofs^2);
                    IA(loc) = set_Rows(idx,nDofs);
                    JA(loc) = set_Cols(p_idx,nDofs);
                    A(loc) = Aloc(:);
                    offset = offset+nDofs^2;
                   
                    locM = offsetM + (1:nDofs^2);
                    IM(locM) = set_Rows(idx,nDofs);
                    JM(locM) = set_Cols(idx,nDofs);
                    M(locM) = Mloc(:);
                    offsetM = offsetM+nDofs^2;
                    
                else
                    % we use exact integration
                    
                    % Preallocate memory
                    Aloc = zeros(nDofs,nDofs);
                    Mloc = Aloc;
                    
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
                        weights = QuadRule.w*Tdet_BK;

                        %evaluate basisfunction of element i at evaluationpoints
                        stevalpoints = (x - ones(nPts,1)*a1)*inv_BK;
                        N = Shap(stevalpoints);

                        % evaluate basisfunction element p_Element in XtauMesh at evaluation point
                        point_hat=(x - ones(nPts,1)*p_bK)*p_inv_BK;
                        p_N =Shap(point_hat);

                        % Add contributions to stiffness matrix
                        for iloc = 1:nDofs
                            for jloc = 1:nDofs
                                Aloc(iloc,jloc) = Aloc(iloc,jloc) + sum(weights.*sum(N(:,iloc).*p_N(:,jloc),2));
                                Mloc(iloc,jloc) = Mloc(iloc,jloc) + sum(weights.*sum(N(:,iloc).*N(:,jloc),2));
                            end
                        end
                    end  % loop over triangles of polygon
                                    
                    % store entries
                    loc = offset + (1:nDofs^2);
                    IA(loc) = set_Rows(idx,nDofs);
                    JA(loc) = set_Cols(p_idx,nDofs);
                    A(loc) = Aloc(:);
                    offset = offset+nDofs^2;
                   
                    locM = offsetM + (1:nDofs^2);
                    IM(locM) = set_Rows(idx,nDofs);
                    JM(locM) = set_Cols(idx,nDofs);
                    M(locM) = Mloc(:);
                    offsetM = offsetM+nDofs^2;
                    
                end % choose quadrature on polygons
            end % if polygon area larger than 0
        end % loop over polygons of element i  
   end % if polygons are not empty 
end % loop over elements of Mesh

% Assign output arguments
NumDOFS =  nElements*nDofs;
              
if(nargout > 2)
    varargout{1} = IA;
    varargout{2} = JA;
    varargout{3} = A;
    varargout{1} = IM;
    varargout{2} = JM;
    varargout{3} = M;
else
    varargout{1} = sparse(IA(1:loc(1)-1),JA(1:loc(1)-1),A(1:loc(1)-1),NumDOFS,NumDOFS);
    varargout{2} = sparse(IM(1:loc(1)-1),JM(1:loc(1)-1),M(1:loc(1)-1),NumDOFS,NumDOFS);
end