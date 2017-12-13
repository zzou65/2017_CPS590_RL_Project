function intersec = aff_elems2(Mesh, defMesh, pbV)
% input:
% Mesh:     original Mesh containing nElemso elements
% defMesh:  deformed Mesh
% pbV:      table containing information about the moved vertices; column 1-2
%           determine the coordinates of every vertex in the deformed mesh;
%           column 3 determines the element number of the original mesh the
%           vertex lies in;
%
% output:
% intersec: 1 by nElems field of structs; struct number i consists of:
%       - Elems: 1 by nElems matrix containing the numbers of the deformed mesh with non
%               0 intersection with element i of the original mesh
%       - nElems: integer number specifying the number of deformed triangles
%               intersecting triangle number i from the original mesh

% algorithm: walk along every edge in the deformed mesh and add the
%   triangles on both sides of the edge to the elements the edge intersects

%   Copyright 2009 Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland


% initialize the constants
%nCoord = size(Mesh.Coordinates,1);
nElemso = size(Mesh.Elements,1);
nEdges = size(Mesh.Edges,1);
tol = eps*10^6;

% preallocate the struct array intersec
intersec(nElemso) = struct('Elems',[],'nElems',0);

% loop over all vertices of the deformed mesh
for i = 1:nEdges
    e1 = defMesh.Edges(i,1);
    e2 = defMesh.Edges(i,2);
    elem1 = pbV(e1,3); % element e1 lies in
    elem2 = pbV(e2,3); % element e2 lies in
    e1coord = defMesh.Coordinates(e1,:);
    e2coord = defMesh.Coordinates(e2,:);
    v1 = Mesh.Coordinates(Mesh.Elements(elem1,1),:); %
    v2 = Mesh.Coordinates(Mesh.Elements(elem1,2),:); % compute vertex coordinates of elem1
    v3 = Mesh.Coordinates(Mesh.Elements(elem1,3),:); %

    loc = defMesh.Edge2Elem(i,:) ~= 0;
    adjelems = defMesh.Edge2Elem(i,loc); % determine the elements sharing the vertex with number loc

    % compute barycentric coordinates of e1 and e2 wrt elem1
    matr = [v1-v3; v2-v3]';
    e1b = matr\(e1coord -v3)';
    e1b = [e1b(1) e1b(2) 1-sum(e1b)];

    e2b = matr\(e2coord -v3)';
    e2b = [e2b(1) e2b(2) 1-sum(e2b)];

    % check if e1 is on the boundary of the triangle (elem1 needs update!)
    if abs(max(e1b) - 1) < tol % e1 is on vertex
        nr = find(abs(e1b - 1) < tol); % find number of vertex e1 lies on
        commvert = Mesh.Elements(elem1,nr);
        adjelements = Mesh.AdjElements(Mesh.Elements(elem1,nr),:); % compute the numbers of all adjacent elements
        nadjelements = Mesh.nAdjElements(Mesh.Elements(elem1,nr),:);
        for j=1:nadjelements
            elem1 = adjelements(j);
            % compute barycentric coordinates of e2 wrt nextelem
            v1 = Mesh.Coordinates(Mesh.Elements(elem1,1),:); %
            v2 = Mesh.Coordinates(Mesh.Elements(elem1,2),:); % compute vertex coordinates of elem1
            v3 = Mesh.Coordinates(Mesh.Elements(elem1,3),:); %
            matr = [v1-v3; v2-v3]';
            e2b = matr\(e2coord -v3)';
            e2b = [e2b(1) e2b(2) 1-sum(e2b)];
            edgenr2 = find(Mesh.Elements(elem1,:)==commvert);
            % barycentric coordinates of e1 wrt nextelem
            e1b = zeros(1,3);
            e1b(edgenr2) = 1;
            t = -e1b(edgenr2)/(e2b(edgenr2)-e1b(edgenr2));
            cb = e1b + t*(e2b - e1b);
            if t > tol && min(cb) >= - tol
                break
            end
        end
    elseif abs(min(e1b)) < tol % e1 is on edge (min(e1b) == 0)
        loc = e1b == 0;
        if e2b(loc) < -tol
            elem1 = Mesh.Neigh(elem1,loc);
        end
    end

    if min(e2b)>=-tol % check if e2 lies already in triangle elem1
        elem2 = elem1;
        vec = abs(e1b) + abs(e2b);
        nr = find(abs(vec) < tol, 1);
        if ~isempty(nr) % check if e1 and e2 lie on the same edge
            neighbor = Mesh.Neigh(elem1,nr); % the neighboring element has to be added also
            %if neighbor ~= -1
            if neighbor > -1  % HH
                intersec(neighbor).Elems = [intersec(neighbor).Elems setdiff(adjelems,intersec(neighbor).Elems)];
            end
        end
    end

    count = 0;

    nextelem = elem1;
    while 1
        intersec(nextelem).Elems = [intersec(nextelem).Elems setdiff(adjelems,intersec(nextelem).Elems)]; % add adjelems to nextelem
        count = count+1;
        if nextelem == elem2
            break; % end of vertex reached
        end
        t1 = -e1b(1)/(e2b(1) - e1b(1));
        c1b = e1b + t1*(e2b - e1b); % intersect e1,e2 with v2,v3
        t2 = -e1b(2)/(e2b(2) - e1b(2));
        c2b = e1b + t2*(e2b - e1b); % intersect e1,e2 with v1,v3
        t3 = -e1b(3)/(e2b(3) - e1b(3));
        c3b = e1b + t3*(e2b - e1b); % intersect e1,e2 with v1,v2

        vertflag = 0;
        if t1 > tol && min(c1b) >=-tol && abs(max(c1b) - 1) < tol % c1b is vertex
            vertflag = 1;
            nr = find(abs(c1b-1) < tol); % find out which local vertex c1b is
            e1coord = c1b(1)*v1 + c1b(2)*v2 + c1b(3)*v3; % update e1coord
        elseif t2 > tol && min(c2b) >=-tol && abs(max(c2b) - 1) < tol
            vertflag = 1;
            nr = find(abs(c2b-1) < tol);
            e1coord = c2b(1)*v1 + c2b(2)*v2 + c2b(3)*v3;
        elseif t3 > tol && min(c3b) >=-tol && abs(max(c3b) - 1) < tol
            vertflag = 1;
            nr = find(abs(c3b-1) < tol);
            e1coord = c3b(1)*v1 + c3b(2)*v2 + c3b(3)*v3;
        elseif (t1 > tol && min(c1b) >= -tol) % intersection point does not lie on a vertex but on an edge
            nextelem = Mesh.Neigh(nextelem,1); % find element that shares the edge edgenr
            e1coord = c1b(1)*v1 + c1b(2)*v2 + c1b(3)*v3; % update e1coord
        elseif (t2 > tol && min(c2b) >= -tol)
            nextelem = Mesh.Neigh(nextelem,2); % find element that shares the edge edgenr
            e1coord = c2b(1)*v1 + c2b(2)*v2 + c2b(3)*v3;
        elseif (t3 > tol && min(c3b) >= -tol)
            nextelem = Mesh.Neigh(nextelem,3); % find element that shares the edge edgenr
            e1coord = c3b(1)*v1 + c3b(2)*v2 + c3b(3)*v3;
        end

        if vertflag % search for element, where the edge continues after leaving nextelem
            commvert = Mesh.Elements(nextelem,nr); % determine the vertex e1coord lies on
            adjelements = Mesh.AdjElements(Mesh.Elements(nextelem,nr),:);
            nadjelements = Mesh.nAdjElements(Mesh.Elements(nextelem,nr),:);
            for j=1:nadjelements
                nextelem = adjelements(j);
                % compute barycentric coordinates of e2 wrt nextelem
                v1 = Mesh.Coordinates(Mesh.Elements(nextelem,1),:); %
                v2 = Mesh.Coordinates(Mesh.Elements(nextelem,2),:); % compute vertex coordinates of nextelem
                v3 = Mesh.Coordinates(Mesh.Elements(nextelem,3),:); %
                matr = [v1-v3; v2-v3]';
                e2b = matr\(e2coord -v3)';
                e2b = [e2b(1) e2b(2) 1-sum(e2b)];
                e1b = matr\(e1coord -v3)';
                e1b = [e1b(1) e1b(2) 1-sum(e1b)];
                edgenr2 = find(Mesh.Elements(nextelem,:)==commvert);
                t = -e1b(edgenr2)/(e2b(edgenr2)-e1b(edgenr2));
                cb = e1b + t*(e2b - e1b); % intersection point with edgenr2
                if t > 0 && min(cb) >=0 % check if the element the edge continues in is nextelem
                    vec = abs(cb) + abs(e1b);
                    nr = find(abs(vec) < tol, 1);
                    if ~isempty(nr) % check if e1b and cb lie on the same edge
                        neighbor = Mesh.Neigh(nextelem,nr);
                        %if neighbor ~= -1 % if neighboring element exists add it also to the intersection
                        if neighbor > -1 % if neighboring element exists add it also to the intersection
                            intersec(neighbor).Elems = [intersec(neighbor).Elems setdiff(adjelems,intersec(neighbor).Elems)];
                        end
                    end
                    break;
                end
            end
        end

        % compute barycentric coordinates of e1 and e2 wrt nextelem
        v1 = Mesh.Coordinates(Mesh.Elements(nextelem,1),:); %
        v2 = Mesh.Coordinates(Mesh.Elements(nextelem,2),:); % compute vertex coordinates of nextelem
        v3 = Mesh.Coordinates(Mesh.Elements(nextelem,3),:); %
        matr = [v1-v3; v2-v3]';
        e1b = matr\(e1coord -v3)';
        e1b = [e1b(1) e1b(2) 1-sum(e1b)];
        e2b = matr\(e2coord -v3)';
        e2b = [e2b(1) e2b(2) 1-sum(e2b)];

        if min(e2b) >= -tol % check if e2b lies in the same element as e1b
            elem2 = nextelem;
        end
    end
end

for i = 1:nElemso % add the number of elements that intersect to every struct
    intersec(i).nElems = length(intersec(i).Elems);
end
