function varargout= assemMat_SemiLagQuad_W1F_patches2ndQuad(Mesh, defMesh,intersec,QuadRule,varargin)
% 
% assemble the stiffness matrix for the semi-lagrange method using one
% point numerical integration on every intersection of triangles of the
% original mesh and the deformed mesh.
%
%   Copyright 2008-2009 Holger Heumann, Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

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
for i = 1:nElements
  polygons = patches(Mesh,defMesh,i,intersec(i));
  for j = 1:intersec(i).nElems  
    
    % compute points of the intersection
    p_Element = polygons(j).Elem;
    polygon = polygons(j).Polygon;
    [polygon,idx,jdx] = unique(polygon,'rows');
    
    if size(polygon,1)>2
        Aloc=zeros(6,6);
        
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
        TK = transpose(inv_BK);

        TRI = delaunay(polygon(:,1),polygon(:,2));
        nTRI = size(TRI,1);
        for k =1:nTRI

            % Compute element mapping
            TbK = polygon(TRI(k,1),:);
            TBK = [polygon(TRI(k,2),:)-TbK; ...
                polygon(TRI(k,3),:)-TbK];
            Tdet_BK = abs(det(TBK));

            x = QuadRule.x*TBK+ones(nPts,1)*TbK;

            % transfer the evalutation point to the standard triangle
            stevalpoints = (x - ones(nPts,1)*a1)*inv_BK;

            % Extract global edge numbers
            eidx = [Mesh.Vert2Edge(Mesh.Elements(i,2),Mesh.Elements(i,3)) ...
                Mesh.Vert2Edge(Mesh.Elements(i,3),Mesh.Elements(i,1)) ...
                Mesh.Vert2Edge(Mesh.Elements(i,1),Mesh.Elements(i,2))];

            N = shap_W1F2nd(stevalpoints);
            N(:,[1 2]) = N(:,[1 2])*TK;
            N(:,[3 4]) = N(:,[3 4])*TK;
            N(:,[5 6]) = N(:,[5 6])*TK;
            N(:,[7 8]) = N(:,[7 8])*TK;
            N(:,[9 10]) = N(:,[9 10])*TK;
            N(:,[11 12]) = N(:,[11 12])*TK;

            % Vertices
            p_vid = Mesh.Elements(p_Element,:);
            p_a1 = defMesh.Coordinates(p_vid(1),:);
            p_a2 = defMesh.Coordinates(p_vid(2),:);
            p_a3 = defMesh.Coordinates(p_vid(3),:);

            p_eidx= [Mesh.Vert2Edge(Mesh.Elements(p_Element,2),Mesh.Elements(p_Element,3)) ...
                Mesh.Vert2Edge(Mesh.Elements(p_Element,3),Mesh.Elements(p_Element,1)) ...
                Mesh.Vert2Edge(Mesh.Elements(p_Element,1),Mesh.Elements(p_Element,2))];

            % Determine the orientation
            if(Mesh.Edges(p_eidx(1),1)==p_vid(2)),  p_p1 = 1;  else    p_p1 = -1;  end
            if(Mesh.Edges(p_eidx(2),1)==p_vid(3)),  p_p2 = 1;  else    p_p2 = -1;  end
            if(Mesh.Edges(p_eidx(3),1)==p_vid(1)),  p_p3 = 1;  else    p_p3 = -1;  end

            % Compute element mapping
            p_bK= p_a1;
            p_BK = [p_a2-p_bK; ...
                p_a3-p_bK];
            p_det_BK = abs(det(p_BK));
            p_inv_BK = inv(p_BK);
            p_TK = transpose(p_inv_BK);

            point_hat=(x - ones(nPts,1)*p_bK)*p_inv_BK;
            p_N =shap_W1F2nd(point_hat);

            p_N(:,[1 2]) = p_N(:,[1 2])*p_TK;
            p_N(:,[3 4]) = p_N(:,[3 4])*p_TK;
            p_N(:,[5 6]) = p_N(:,[5 6])*p_TK;
            p_N(:,[7 8]) = p_N(:,[7 8])*p_TK;
            p_N(:,[9 10]) = p_N(:,[9 10])*p_TK;
            p_N(:,[11 12]) = p_N(:,[11 12])*p_TK;

            % Determine the orientation
            if(Mesh.Edges(eidx(1),1)==vid(2)),  p1 = 1;  else    p1 = -1;  end
            if(Mesh.Edges(eidx(2),1)==vid(3)),  p2 = 1;  else    p2 = -1;  end
            if(Mesh.Edges(eidx(3),1)==vid(1)),  p3 = 1;  else    p3 = -1;  end

            w = QuadRule.w'*Tdet_BK;
            % Add contributions to stiffness matrix
            Aloc(1:3,1:3) = Aloc(1:3,1:3)+...
                [p1*sum(w*sum(N(:,[1 2]).*p_N(:,[1 2]),2))*p_p1 p1*sum(w*sum(N(:,[1 2]).*p_N(:,[3 4]),2))*p_p2 p1*sum(w*sum(N(:,[1 2]).*p_N(:,[5 6]),2))*p_p3; ...
                p2*sum(w*sum(N(:,[3 4]).*p_N(:,[1 2]),2))*p_p1 p2*sum(w*sum(N(:,[3 4]).*p_N(:,[3 4]),2))*p_p2 p2*sum(w*sum(N(:,[3 4]).*p_N(:,[5 6]),2))*p_p3; ...
                p3*sum(w*sum(N(:,[5 6]).*p_N(:,[1 2]),2))*p_p1 p3*sum(w*sum(N(:,[5 6]).*p_N(:,[3 4]),2))*p_p2 p3*sum(w*sum(N(:,[5 6]).*p_N(:,[5 6]),2))*p_p3];

            Aloc(1:3,4:6) = Aloc(1:3,4:6)+...
                [p1*sum(w*sum(N(:,[1 2]).*p_N(:,[7 8]),2)) p1*sum(w*sum(N(:,[1 2]).*p_N(:,[9 10]),2)) p1*sum(w*sum(N(:,[1 2]).*p_N(:,[11 12]),2)); ...
                p2*sum(w*sum(N(:,[3 4]).*p_N(:,[7 8]),2)) p2*sum(w*sum(N(:,[3 4]).*p_N(:,[9 10]),2)) p2*sum(w*sum(N(:,[3 4]).*p_N(:,[11 12]),2)); ...
                p3*sum(w*sum(N(:,[5 6]).*p_N(:,[7 8]),2)) p3*sum(w*sum(N(:,[5 6]).*p_N(:,[9 10]),2)) p3*sum(w*sum(N(:,[5 6]).*p_N(:,[11 12]),2))];

            Aloc(4:6,1:3)= Aloc(4:6,1:3)+...
                [sum(w*sum(N(:,[7 8]).*p_N(:,[1 2]),2))*p_p1     sum(w*sum(N(:,[7 8]).*p_N(:,[3 4]),2))*p_p2      sum(w*sum(N(:,[7 8]).*p_N(:,[5 6]),2))*p_p3; ...
                sum(w*sum(N(:,[9 10]).*p_N(:,[1 2]),2))*p_p1   sum(w*sum(N(:,[9 10]).*p_N(:,[3 4]),2))*p_p2    sum(w*sum(N(:,[9 10]).*p_N(:,[5 6]),2))*p_p3; ...
                sum(w*sum(N(:,[11 12]).*p_N(:,[1 2]),2))*p_p1 sum(w*sum(N(:,[11 12]).*p_N(:,[3 4]),2))*p_p2  sum(w*sum(N(:,[11 12]).*p_N(:,[5 6]),2))*p_p3];

            Aloc(4:6, 4:6) = Aloc(4:6, 4:6)+...
                [sum(w*sum(N(:,[7 8]).*p_N(:,[7 8]),2))    sum(w*sum(N(:,[7 8]).*p_N(:,[9 10]),2))      sum(w*sum(N(:,[7 8]).*p_N(:,[11 12]),2)); ...
                sum(w*sum(N(:,[9 10]).*p_N(:,[7 8]),2))   sum(w*sum(N(:,[9 10]).*p_N(:,[9 10]),2))   sum(w*sum(N(:,[9 10]).*p_N(:,[11 12]),2)); ...
                sum(w*sum(N(:,[11 12]).*p_N(:,[7 8]),2)) sum(w*sum(N(:,[11 12]).*p_N(:,[9 10]),2)) sum(w*sum(N(:,[11 12]).*p_N(:,[11 12]),2))];

        end  % loop over triangles of polygon
        I(loc) = set_Rows([eidx, nEdges+eidx],6);
        J(loc) = set_Cols([p_eidx, nEdges+p_eidx],6);
        A(loc) = Aloc(:);
        loc = loc+36;
    end  % if size(polygon,1)>2 else do nothing
  end % loop over intersetions
end % loop over elements

% Assign output arguments
if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
else
    varargout{1} = sparse(I(1:loc-36),J(1:loc-36),A(1:loc-36),2*nEdges,2*nEdges);
end
