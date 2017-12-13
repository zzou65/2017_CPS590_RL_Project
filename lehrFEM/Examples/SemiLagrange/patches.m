function polygons = patches(Mesh, defMesh, elem, intersec)
% computes the boundary of the intersection of an element of the original
% mesh with all the intersecting triangles of the deformed mesh

% input:
% Mesh:     original mesh
% defMesh:  deformed mesh
% elem:     element number specifying the element the polygon should be
%           computed
% intersec: struct containing the fields
%           - Elems: 1 by nElems matrix containing the elements with non zero intersection
%           - nElems:number specifying the amount of triangles intersecting
%               with elem 

%   Copyright 2009 Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% determine the tolerance
tol = eps*10^6;

% preallocate space
polygons(intersec.nElems) = struct('Elem',0,'Polygon',[]);

% determine the vertices of elem
v1 = Mesh.Coordinates(Mesh.Elements(elem,1),:);
v2 = Mesh.Coordinates(Mesh.Elements(elem,2),:);
v3 = Mesh.Coordinates(Mesh.Elements(elem,3),:);
v(1,:) = v1;
v(2,:) = v2;
v(3,:) = v3;

% if only one triangle of the deformed mesh intersects the original element
% elem is covered fully by it.

if intersec.nElems == 1
    polygons(1).Polygon = [v1;v2;v3];
    polygons(1).Elem = intersec.Elems(1);
    return
end

% compute the matrix for computing barycentric coordinates
matr = [v1-v3; v2-v3]';

for i=1:intersec.nElems
    interelem = intersec.Elems(i);
    % compute coordinates of the vertices of the deformed triangle
    w1 = defMesh.Coordinates(defMesh.Elements(interelem,1),:);
    w2 = defMesh.Coordinates(defMesh.Elements(interelem,2),:);
    w3 = defMesh.Coordinates(defMesh.Elements(interelem,3),:);
    
  
    
    % find the coordinates that lie inside the original triangle
    wb1 = matr\(w1-v3)';
    wb2 = matr\(w2-v3)';
    wb3 = matr\(w3-v3)';
    wb1 = [wb1' 1-sum(wb1)];
    wb2 = [wb2' 1-sum(wb2)];
    wb3 = [wb3' 1-sum(wb3)];
    
    % all coordinates are inside the triangle
    if (min(wb1) >= -tol && min(wb2) >= -tol && min(wb3) >= -tol)
        polygons(i).Polygon = [w1;w2;w3];
        polygons(i).Elem = intersec.Elems(i);
    end
    
    % two coordinates are inside the triangle
    twoflag = 0; % 1 if two vertices of deformed element are in elem, 0 otherwise
    if (min(wb1) >= -eps && min(wb2) >= -tol)
        twoflag = 1;
    elseif (min(wb2) >= -tol && min(wb3) >= -tol)
        twoflag = 1;
        temp = wb1;  %
        wb1 = wb3;   % renumbering of vertices wbi to have wb1 and wb2 inside
        wb3 = temp;  %
        
        temp = w1;  %
        w1 = w3;    % renumbering of vertices wbi to have wb1 and wb2 inside
        w3 = temp;  %
    elseif (min(wb1) >= -tol && min(wb3) >= -tol)
        twoflag = 1;
        temp = wb2;  %
        wb2 = wb3;   % renumbering of vertices wbi to have wb1 and wb2 inside
        wb3 = temp;  %
        
        temp = w2;  %
        w2 = w3;    % renumbering of vertices wbi to have wb1 and wb2 inside
        w3 = temp;  %
    end
    
    if twoflag % two vertices lie inside the triangle
        % compute the intersection points of w1,w3 with edges of elem
        t1 = -wb1(1)/(wb3(1)-wb1(1));
        t2 = -wb1(2)/(wb3(2)-wb1(2));
        t3 = -wb1(3)/(wb3(3)-wb1(3));
        y(1,:) = wb1 + t1*(wb3-wb1); % intersection with v2,v3
        y(2,:) = wb1 + t2*(wb3-wb1); % intersection with v1,v3
        y(3,:) = wb1 + t3*(wb3-wb1); % intersection with v1,v2
        inters1 = 0;
        if (t1 > 0 && min(y(1,[2 3])) >= -tol) % intersection with v2,v3
            inters1 = 1;
        elseif (t2 > 0 && min(y(2,[1 3])) >= -tol) % intersection with v1,v3
            inters1 = 2;
        elseif (t3 > 0 && min(y(3,[1 2])) >= -tol) % intersection with v1,v2
            inters1 = 3;
        end
        
        % compute the intersection points of w2,w3 with edges of elem
        t1 = -wb2(1)/(wb3(1)-wb2(1));
        t2 = -wb2(2)/(wb3(2)-wb2(2));
        t3 = -wb2(3)/(wb3(3)-wb2(3));
        z(1,:) = wb2 + t1*(wb3-wb2); % intersection with v2,v3
        z(2,:) = wb2 + t2*(wb3-wb2); % intersection with v1,v3
        z(3,:) = wb2 + t3*(wb3-wb2); % intersection with v1,v2
        inters2=0;
        if (t1 > 0 && min(z(1,[2 3])) >= -tol) % intersection with v2,v3
            inters2 = 1;
        elseif (t2 > 0 && min(z(2,[1 3])) >= -tol) % intersection with v1,v3
            inters2 = 2;
        elseif (t3 > 0 && min(z(3,[1 2])) >= -tol) % intersection with v1,v2
            inters2 = 3;
        end
        
        % assemble the polygons from the intersection points
        vec = abs(wb1) + abs(wb2);
        nr = find(vec < tol, 1);
        if ~isempty(nr) % w1 and w2 lie on the same edge
            if wb3(nr) < 0 % check on which side wb3 lies
                polygons(i).Polygon = [w1;w2];
                polygons(i).Elem = intersec.Elems(i);
            elseif min(wb3) >= -tol % w3 lies inside the triangle
                polygons(i).Polygon = [w1;w2;w3];
                polygons(i).Elem = intersec.Elems(i);
            else
                if inters2 % w2,w3 intersects with one of the vertices
                    zz = z(inters2,1)*v1 + z(inters2,2)*v2 + z(inters2,3)*v3; 
                    polygons(i).Polygon = [w1;w2;zz];
                    polygons(i).Elem = intersec.Elems(i);
                elseif inters1 % w1,w3 intersects with one of the vertices
                    yy = y(inters1,1)*v1 + y(inters1,2)*v2 + y(inters1,3)*v3;
                    polygons(i).Polygon = [w1;w2;yy];
                    polygons(i).Elem = intersec.Elems(i);
                else % no intersections, w1,w2,w3 covers the whole triangle
                    polygons(i).Polygon = [v1;v2;v3];
                    polygons(i).Elem = intersec.Elems(i);
                end
            end
        elseif ~inters2 % w2,w3 does not intersect with any edge
            yy = y(inters1,1)*v1 + y(inters1,2)*v2 + y(inters1,3)*v3;
            polygons(i).Polygon = [w1;w2;yy];
            polygons(i).Elem = intersec.Elems(i);
        elseif ~inters1 % w1,w3 does not intersect with any edge
            zz = z(inters2,1)*v1 + z(inters2,2)*v2 + z(inters2,3)*v3;
            polygons(i).Polygon = [w1;w2;zz];
            polygons(i).Elem = intersec.Elems(i);
        elseif inters1 == inters2 % w1,w3 and w2,w3 intersect the same edge
            zz = z(inters1,1)*v1 + z(inters1,2)*v2 + z(inters1,3)*v3;
            yy = y(inters1,1)*v1 + y(inters1,2)*v2 + y(inters1,3)*v3;
            polygons(i).Polygon = [w1;w2;zz;yy];
            polygons(i).Elem = intersec.Elems(i);
        else % w1,w3 and w2,w3 intersect different edges
            midvert = setdiff([1 2 3], [inters1 inters2]); % determine the middle node
            midcoord = v(midvert,:);
            zz = z(inters2,1)*v1 + z(inters2,2)*v2 + z(inters2,3)*v3;
            yy = y(inters1,1)*v1 + y(inters1,2)*v2 + y(inters1,3)*v3;
            polygons(i).Polygon = [w1;w2;zz;midcoord;yy];
            polygons(i).Elem = intersec.Elems(i);
        end
    end
    
    
    % one coordinate is inside the triangle
    oneflag = 0;
    if min(wb1) >= -tol
        oneflag = 1;
    elseif min(wb2) >= -tol % renumber vertices to have w1 inside
        temp = wb1;
        wb1 = wb2;
        wb2 = temp;
        
        temp2 = w1;
        w1 = w2;
        w2 = temp2;
        oneflag = 1;
    elseif min(wb3) >= -tol % renumber vertices to have w1 inside
        temp = wb1;
        wb1 = wb3;
        wb3 = temp;
        
        temp2 = w1;
        w1 = w3;
        w3 = temp2;
        oneflag = 1;
    end
    
    if oneflag && ~twoflag % one vertex lies inside
        % compute intersection of w1,w2 with the three edges of elem
        t1 = -wb1(1)/(wb2(1)-wb1(1));
        t2 = -wb1(2)/(wb2(2)-wb1(2));
        t3 = -wb1(3)/(wb2(3)-wb1(3));
        y(1,:) = wb1 + t1*(wb2-wb1); % intersection with v2,v3
        y(2,:) = wb1 + t2*(wb2-wb1); % intersection with v1,v3
        y(3,:) = wb1 + t3*(wb2-wb1); % intersection with v1,v2
        inters1 = 0;
        if (t1 > tol && min(y(1,[2 3])) >= -tol) % intersection with v2,v3
            inters1 = 1;
        elseif (t2 > tol && min(y(2,[1 3])) >= -tol) % intersection with v1,v3
            inters1 = 2;
        elseif (t3 > tol && min(y(3,[1 2])) >= -tol) % intersection with v1,v2
            inters1 = 3;
        end
        
        % compute intersection of w1,w3 with the three edges of elem
        s1 = -wb1(1)/(wb3(1)-wb1(1));
        s2 = -wb1(2)/(wb3(2)-wb1(2));
        s3 = -wb1(3)/(wb3(3)-wb1(3));
        z(1,:) = wb1 + s1*(wb3-wb1); % intersection with v2,v3
        z(2,:) = wb1 + s2*(wb3-wb1); % intersection with v1,v3
        z(3,:) = wb1 + s3*(wb3-wb1); % intersection with v1,v2
        inters2 = 0;
        if (s1 > tol && min(z(1,[2 3])) >= -tol) % intersection with v2,v3
            inters2 = 1;
        elseif (s2 > tol && min(z(2,[1 3])) >= -tol) % intersection with v1,v3
            inters2 = 2;
        elseif (s3 > tol && min(z(3,[1 2])) >= -tol) % intersection with v1,v2
            inters2 = 3;
        end
        
        % assemble the polygons and process the different cases occuring
        if ~inters1 && ~inters2 % no intersections (w1 lies on vertex of original triangle)
            if norm(w1-v1) < tol
                
            elseif norm(w1-v2) < tol % renumber the original vertices
                temp = v2;
                v2 = v1;
                v1 = temp;
                temp = wb1(2);
                wb1(2)=wb1(1);
                wb1(1)=temp;
                temp = wb2(2);
                wb2(2)=wb2(1);
                wb2(1)=temp;
                temp = wb3(2);
                wb3(2)=wb3(1);
                wb3(1)=temp;
            elseif norm(w1-v3) < tol % renumber the original vertices
                temp = v3;
                v3 = v1;
                v1 = temp;
                temp = wb1(3);
                wb1(3)=wb1(1);
                wb1(1)=temp;
                temp = wb2(3);
                wb2(3)=wb2(1);
                wb2(1)=temp;
                temp = wb3(3);
                wb3(3)=wb3(1);
                wb3(1)=temp;
            end
            
            % compute intersection of w2,w3 with the three edges of elem
            r1 = -wb2(1)/(wb3(1)-wb2(1));
            r2 = -wb2(2)/(wb3(2)-wb2(2));
            r3 = -wb2(3)/(wb3(3)-wb2(3));
            x(1,:) = wb2 + r1*(wb3-wb2); % intersection with v2,v3
            x(2,:) = wb2 + r2*(wb3-wb2); % intersection with v1,v3
            x(3,:) = wb2 + r3*(wb3-wb2); % intersection with v1,v2
            
            xx1 = x(1,1)*v1 + x(1,2)*v2 + x(1,3)*v3;
            xx2 = x(2,1)*v1 + x(2,2)*v2 + x(2,3)*v3;
            xx3 = x(3,1)*v1 + x(3,2)*v2 + x(3,3)*v3;
           
            
            if r2 > 0 && min(x(2,:)) >= -tol && r3 > 0 && min(x(3,:)) >= -tol
                polygons(i).Polygon = [w1;xx2;xx3];
                polygons(i).Elem = intersec.Elems(i);
            elseif r2 > 0 && min(x(2,:)) >= -tol && r1 > 0 && min(x(1,:)) >= -tol
                polygons(i).Polygon = [w1;xx2;xx1;v2];
                polygons(i).Elem = intersec.Elems(i);
            elseif r3 > 0 && min(x(3,:)) >= -tol && r1 > 0 && min(x(1,:)) >= -tol
                polygons(i).Polygon = [w1;xx3;xx1;v3];
                polygons(i).Elem = intersec.Elems(i);
            else
                polygons(i).Polygon = w1;
                polygons(i).Elem = intersec.Elems(i);
            end
                        
        elseif ~inters1 % no intersection of w1,w2 with the vertices of the element
            loc = setdiff([1 2 3], inters2);
            vv1 = v(loc(1),:);
            vv2 = v(loc(2),:);
            coord3 = setdiff([1 2 3],loc);
            zz = z(inters2,1)*v1 + z(inters2,2)*v2 + z(inters2,3)*v3;
            
            % intersect vv1,vv2 with w2,w3
            t = -wb2(coord3)/(wb3(coord3)-wb2(coord3));
            pointb = wb2 + t*(wb3-wb2);
            point = pointb(1)*v1 + pointb(2)*v2 + pointb(3)*v3;
            if norm(w2 - vv1) <= norm(w2 - vv2) % check on which side w2 lies
                if t>0 && min(pointb) >= -tol % check if point is inside the triangle
                    t = -wb2(loc(2))/(wb3(loc(2))-wb2(loc(2)));
                    point2b = wb2 + t*(wb3-wb2);
                    point2 = point2b(1)*v1 + point2b(2)*v2 + point2b(3)*v3;
                    polygons(i).Polygon = [w1;zz;point;point2];
                    polygons(i).Elem = intersec.Elems(i);
                else
                    polygons(i).Polygon = [w1;zz;vv1];
                    polygons(i).Elem = intersec.Elems(i);
                end
            else
                if t>0 && min(pointb) >= -tol  % check if point is inside the triangle
                    t = -wb2(loc(1))/(wb3(loc(1))-wb2(loc(1)));
                    point2b = wb2 + t*(wb3-wb2);
                    point2 = point2b(1)*v1 + point2b(2)*v2 + point2b(3)*v3;
                    polygons(i).Polygon = [w1;zz;point;point2];
                    polygons(i).Elem = intersec.Elems(i);
                else
                    polygons(i).Polygon = [w1;zz;vv2];
                    polygons(i).Elem = intersec.Elems(i);
                end
            end
            
        elseif ~inters2  % no intersection of w1,w3 with the vertices of the element
            loc = setdiff([1 2 3], inters1);
            vv1 = v(loc(1),:);
            vv2 = v(loc(2),:);
            coord3 = setdiff([1 2 3],loc);
            vv3 = v(coord3,:);
            yy = y(inters1,1)*v1 + y(inters1,2)*v2 + y(inters1,3)*v3;
            
            % intersect vv1,vv2 with w2,w3
            t = -wb2(coord3)/(wb3(coord3)-wb2(coord3));
            pointb = wb2 + t*(wb3-wb2);
            point = pointb(1)*v1 + pointb(2)*v2 + pointb(3)*v3;
            if norm(w3 - vv1) <= norm(w3 - vv2) % check on which side w3 lies
                if t>0 && min(pointb) >= -tol  % check if point is inside the triangle
                    t = -wb2(loc(2))/(wb3(loc(2))-wb2(loc(2)));
                    point2b = wb2 + t*(wb3-wb2);
                    point2 = point2b(1)*v1 + point2b(2)*v2 + point2b(3)*v3;
                    polygons(i).Polygon = [w1;yy;point;point2];
                    polygons(i).Elem = intersec.Elems(i);
                else
                    polygons(i).Polygon = [w1;yy;vv1];
                    polygons(i).Elem = intersec.Elems(i);
                end
            else
                if t>0 && min(pointb) >= -tol  % check if point is inside the triangle
                    t = -wb2(loc(1))/(wb3(loc(1))-wb2(loc(1)));
                    point2b = wb2 + t*(wb3-wb2);
                    point2 = point2b(1)*v1 + point2b(2)*v2 + point2b(3)*v3;
                    polygons(i).Polygon = [w1;yy;point;point2];
                    polygons(i).Elem = intersec.Elems(i);
                else
                    polygons(i).Polygon = [w1;yy;vv2];
                    polygons(i).Elem = intersec.Elems(i);
                end
            end
            
        elseif inters1 == inters2 % w1,w2 and w1,w3 intersect the same edge
            zz = z(inters2,1)*v1 + z(inters2,2)*v2 + z(inters2,3)*v3;
            yy = y(inters1,1)*v1 + y(inters1,2)*v2 + y(inters1,3)*v3;
            polygons(i).Polygon = [w1;zz;yy];
            polygons(i).Elem = intersec.Elems(i);
        else % w1,w2 and w1,w3 intersect different edges
            % determine intersection point of w2,w3 with edge inters1
            t1 = -wb2(inters1)/(wb3(inters1)-wb2(inters1));
            cb1 = wb2 + t1*(wb3-wb2);
            if min(cb1) >= -tol % check if cb1 lies inside the triangle
                % determine intersection point of w2,w3 with edge inters2
                t2 = -wb2(inters2)/(wb3(inters2)-wb2(inters2));
                cb2 = wb2 + t2*(wb3-wb2);
                c1 = cb1(1)*v1 + cb1(2)*v2 + cb1(3)*v3;
                c2 = cb2(1)*v1 + cb2(2)*v2 + cb2(3)*v3;
                zz = z(inters2,1)*v1 + z(inters2,2)*v2 + z(inters2,3)*v3;
                yy = y(inters1,1)*v1 + y(inters1,2)*v2 + y(inters1,3)*v3;
                polygons(i).Polygon = [w1;yy;c1;c2;zz];
                polygons(i).Elem = intersec.Elems(i);
            else
                midvert = setdiff([1 2 3], [inters1 inters2]); % determine the middle node
                midcoord = v(midvert,:);
                zz = z(inters2,1)*v1 + z(inters2,2)*v2 + z(inters2,3)*v3;
                yy = y(inters1,1)*v1 + y(inters1,2)*v2 + y(inters1,3)*v3;
                polygons(i).Polygon = [w1;zz;midcoord;yy];
                polygons(i).Elem = intersec.Elems(i);
            end
        end
    end
    
    % all coordinates are outside the triangle
    if ~oneflag && ~twoflag
        flag1 = 0; flag2 = 0; flag3 = 0;
        
        % compute intersection of w1,w2 with the three edges of elem
        t1 = -wb1(1)/(wb2(1)-wb1(1));
        t2 = -wb1(2)/(wb2(2)-wb1(2));
        t3 = -wb1(3)/(wb2(3)-wb1(3));
        y(1,:) = wb1 + t1*(wb2-wb1); % intersection with v2,v3
        y(2,:) = wb1 + t2*(wb2-wb1); % intersection with v1,v3
        y(3,:) = wb1 + t3*(wb2-wb1); % intersection with v1,v2
        t = [t1 t2 t3];
        
        loc3 = (t>0) + (t<1);
        loc3 = find(loc3 == 2);
        
        if length(loc3) == 2 && min(min(y(loc3,:))) >= -tol % two intersections
            flag1 = 1;
        elseif length(loc3) == 3 % 3 intersections
            loc3 = [];
            if min(y(1,:)) >= -tol
                loc3 = [loc3 1];
            end
            if min(y(2,:)) >= -tol
                loc3 = [loc3 2];
            end
            if min(y(3,:)) >= -tol
                loc3 = [loc3 3];
            end
            if length(loc3)==2
                flag1=1;
            end
        end
        
        % compute intersection of w1,w3 with the three edges of elem
        s1 = -wb1(1)/(wb3(1)-wb1(1));
        s2 = -wb1(2)/(wb3(2)-wb1(2));
        s3 = -wb1(3)/(wb3(3)-wb1(3));
        z(1,:) = wb1 + s1*(wb3-wb1); % intersection with v2,v3
        z(2,:) = wb1 + s2*(wb3-wb1); % intersection with v1,v3
        z(3,:) = wb1 + s3*(wb3-wb1); % intersection with v1,v2
        s = [s1 s2 s3];
        
        loc2 = (s>0) + (s<1);
        loc2 = find(loc2 == 2);
        
        
        if length(loc2) == 2 && min(min(z(loc2,:))) >= -tol % two intersections
            flag2 = 1;
        elseif length(loc2) == 3 % 3 intersections
            loc2 = [];
            if min(z(1,:)) >= -tol
                loc2 = [loc2 1];
            end
            if min(z(2,:)) >= -tol
                loc2 = [loc2 2];
            end
            if min(z(3,:)) >= -tol
                loc2 = [loc2 3];
            end
            if length(loc2)==2
                flag2=1;
            end
        end
        
        % compute intersection of w2,w3 with the three edges of elem
        r1 = -wb2(1)/(wb3(1)-wb2(1));
        r2 = -wb2(2)/(wb3(2)-wb2(2));
        r3 = -wb2(3)/(wb3(3)-wb2(3));
        x(1,:) = wb2 + r1*(wb3-wb2); % intersection with v2,v3
        x(2,:) = wb2 + r2*(wb3-wb2); % intersection with v1,v3
        x(3,:) = wb2 + r3*(wb3-wb2); % intersection with v1,v2
        r = [r1 r2 r3];
        
        loc1 = (r>0) + (r<1);
        loc1 = find(loc1 == 2);
        
        if length(loc1) >= 2 && min(min(x(loc1,:))) >= -tol % two intersections
            flag3 = 1; % 3 intersections
        elseif length(loc1) == 3
            loc1 = [];
            if min(x(1,:)) >= -tol
                loc1 = [loc1 1];
            end
            if min(x(2,:)) >= -tol
                loc1 = [loc1 2];
            end
            if min(x(3,:)) >= -tol
                loc1 = [loc1 3];
            end
            if length(loc1)==2
                flag3=1;
            end
        end
        
        flag = [flag1,flag2,flag3]; % check which of the vertices intersects
        if sum(flag) == 1 % one edge intersects the original triangle
            nr = find(flag == 1);
            % find the barycentric coordinates of v1,v2,v3 wrt wi
            matr2 = [w1-w3; w2-w3]';
            vb1 = matr2\(v1-w3)';
            vb2 = matr2\(v2-w3)';
            vb3 = matr2\(v3-w3)';
            vvb(1,:) = [vb1' 1-sum(vb1)];
            vvb(2,:) = [vb2' 1-sum(vb2)];
            vvb(3,:) = [vb3' 1-sum(vb3)];
            if nr == 1
                loc = loc3;
                vertnr = setdiff([1 2 3], loc);
                yy1 = y(loc(1),1)*v1 + y(loc(1),2)*v2 + y(loc(1),3)*v3;
                yy2 = y(loc(2),1)*v1 + y(loc(2),2)*v2 + y(loc(2),3)*v3;
                if vvb(vertnr,3) > -tol
                    polygons(i).Polygon = [yy1;yy2;v(vertnr,:)];
                    polygons(i).Elem = intersec.Elems(i);
                else
                    polygons(i).Polygon = [yy1;yy2;v(loc(1),:);v(loc(2),:)];
                    polygons(i).Elem = intersec.Elems(i);
                end
            elseif nr == 2
                loc = loc2;
                vertnr = setdiff([1 2 3], loc);
                zz1 = z(loc(1),1)*v1 + z(loc(1),2)*v2 + z(loc(1),3)*v3;
                zz2 = z(loc(2),1)*v1 + z(loc(2),2)*v2 + z(loc(2),3)*v3;
                if vvb(vertnr,2) > -tol
                    polygons(i).Polygon = [zz1;zz2;v(vertnr,:)];
                    polygons(i).Elem = intersec.Elems(i);
                else
                    polygons(i).Polygon = [zz1;zz2;v(loc(1),:);v(loc(2),:)];
                    polygons(i).Elem = intersec.Elems(i);
                end
            elseif nr == 3
                loc =loc1;
                vertnr = setdiff([1 2 3], loc);
                xx1 = x(loc(1),1)*v1 + x(loc(1),2)*v2 + x(loc(1),3)*v3;
                xx2 = x(loc(2),1)*v1 + x(loc(2),2)*v2 + x(loc(2),3)*v3;
                if vvb(vertnr,1) > -tol
                    polygons(i).Polygon = [xx1;xx2;v(vertnr,:)];
                    polygons(i).Elem = intersec.Elems(i);
                else
                    polygons(i).Polygon = [xx1;xx2;v(loc(1),:);v(loc(2),:)];
                    polygons(i).Elem = intersec.Elems(i);
                end
            end
            
        elseif sum(flag) == 2 % two edges intersect the original triangle
            if flag(1) == 1 && flag(2) == 1
                vertnr = intersect(loc3,loc2);
                if length(vertnr) == 1
                    zz1 = z(vertnr,1)*v1 + z(vertnr,2)*v2 + z(vertnr,3)*v3;
                    zz2 = z(setdiff(loc2,vertnr),1)*v1 + z(setdiff(loc2,vertnr),2)*v2 + z(setdiff(loc2,vertnr),3)*v3;
                    yy1 = y(vertnr,1)*v1 + y(vertnr,2)*v2 + y(vertnr,3)*v3;
                    yy2 = y(setdiff(loc3,vertnr),1)*v1 + y(setdiff(loc3,vertnr),2)*v2 + y(setdiff(loc3,vertnr),3)*v3;
                    polygons(i).Polygon = [yy1;yy2;v(vertnr,:);zz2;zz1];
                    polygons(i).Elem = intersec.Elems(i); 
                else
                    loc2 = [loc2(2) loc2(1)];
                    yy12 = y(loc3,1)*v1 + y(loc3,2)*v2 + y(loc3,3)*v3;
                    zz12 = z(loc2,1)*v1 + z(loc2,2)*v2 + z(loc2,3)*v3;
                    polygons(i).Polygon = [yy12;zz12];
                    polygons(i).Elem = intersec.Elems(i); 
                end
            elseif flag(1) == 1 && flag(3) == 1
                vertnr = intersect(loc3,loc1);
                if length(vertnr) == 1
                    xx1 = x(vertnr,1)*v1 + x(vertnr,2)*v2 + x(vertnr,3)*v3;
                    xx2 = x(setdiff(loc1,vertnr),1)*v1 + x(setdiff(loc1,vertnr),2)*v2 + x(setdiff(loc1,vertnr),3)*v3;
                    yy1 = y(vertnr,1)*v1 + y(vertnr,2)*v2 + y(vertnr,3)*v3;
                    yy2 = y(setdiff(loc3,vertnr),1)*v1 + y(setdiff(loc3,vertnr),2)*v2 + y(setdiff(loc3,vertnr),3)*v3;
                    polygons(i).Polygon = [yy1;yy2;v(vertnr,:);xx2;xx1];
                    polygons(i).Elem = intersec.Elems(i); 
                else
                    loc1 = [loc1(2) loc1(1)];
                    yy12 = y(loc3,1)*v1 + y(loc3,2)*v2 + y(loc3,3)*v3;
                    xx12 = x(loc1,1)*v1 + x(loc1,2)*v2 + x(loc1,3)*v3;
                    polygons(i).Polygon = [yy12;xx12];
                    polygons(i).Elem = intersec.Elems(i); 
                end
            elseif flag(2) == 1 && flag(3) == 1
                vertnr = intersect(loc2,loc1);
                if length(vertnr) == 1
                    zz1 = z(vertnr,1)*v1 + z(vertnr,2)*v2 + z(vertnr,3)*v3;
                    zz2 = z(setdiff(loc2,vertnr),1)*v1 + z(setdiff(loc2,vertnr),2)*v2 + z(setdiff(loc2,vertnr),3)*v3;
                    xx1 = x(vertnr,1)*v1 + x(vertnr,2)*v2 + x(vertnr,3)*v3;
                    xx2 = x(setdiff(loc1,vertnr),1)*v1 + x(setdiff(loc1,vertnr),2)*v2 + x(setdiff(loc1,vertnr),3)*v3;
                    polygons(i).Polygon = [xx1;xx2;v(vertnr,:);zz2;zz1];
                    polygons(i).Elem = intersec.Elems(i); 
                else
                    loc1 = [loc1(2) loc1(1)];
                    zz12 = z(loc2,1)*v1 + z(loc2,2)*v2 + z(loc2,3)*v3;
                    xx12 = x(loc1,1)*v1 + x(loc1,2)*v2 + x(loc1,3)*v3;
                    polygons(i).Polygon = [zz12;xx12];
                    polygons(i).Elem = intersec.Elems(i); 
                end
            end
            
        elseif sum(flag) == 3 % three edges intersects the original triangle
            yy1 = y(loc3(1),1)*v1 + y(loc3(1),2)*v2 + y(loc3(1),3)*v3;
            yy2 = y(loc3(2),1)*v1 + y(loc3(2),2)*v2 + y(loc3(2),3)*v3;
            nr = find(loc2 == loc3(2));
            if ~isempty(nr)
                nextedge = loc2(nr);
                zz1 = z(nextedge,1)*v1 + z(nextedge,2)*v2 + z(nextedge,3)*v3;
                zz2 = z(setdiff(loc2,nextedge),1)*v1 + z(setdiff(loc2,nextedge),2)*v2 + z(setdiff(loc2,nextedge),3)*v3;
                nr = find(loc1 == setdiff(loc2,nextedge));
                nnextedge = loc1(nr);
                xx1 = x(nnextedge,1)*v1 + x(nnextedge,2)*v2 + x(nnextedge,3)*v3;
                xx2 = x(setdiff(loc1,nnextedge),1)*v1 + x(setdiff(loc1,nnextedge),2)*v2 + x(setdiff(loc1,nnextedge),3)*v3;
            else
                temp = yy1;
                yy1 = yy2;
                yy2 = temp;
                nr = find(loc2 == loc3(1));
                nextedge = loc3(nr);
                zz1 = z(nextedge,1)*v1 + z(nextedge,2)*v2 + z(nextedge,3)*v3;
                zz2 = z(setdiff(loc2,nextedge),1)*v1 + z(setdiff(loc2,nextedge),2)*v2 + z(setdiff(loc2,nextedge),3)*v3;
                nr = find(loc1 == setdiff(loc2,nextedge));
                nnextedge = loc1(nr);
                xx1 = x(nnextedge,1)*v1 + x(nnextedge,2)*v2 + x(nnextedge,3)*v3;
                xx2 = x(setdiff(loc1,nnextedge),1)*v1 + x(setdiff(loc1,nnextedge),2)*v2 + x(setdiff(loc1,nnextedge),3)*v3;
            end
            polygons(i).Polygon = [yy1;yy2;zz1;zz2;xx1;xx2];
            polygons(i).Elem = intersec.Elems(i); 
        end
           
    end
end

for i = 1:intersec.nElems % remove double points in the polygons
    nPoints = size(polygons(i).Polygon,1);
    coords = [];
    for j = 0:nPoints-1
        if norm(polygons(i).Polygon(mod(j,nPoints)+1,:)-polygons(i).Polygon(mod(j+1,nPoints)+1,:))==0
            coords = [coords, j+1];
        end
    end
    polygons(i).Polygon(coords,:) = [];
end





