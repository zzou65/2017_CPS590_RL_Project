function Aloc=trace_vertex(vertex,vertex_r,vertex_q,r_id,q_id,gradN,Element,velo,Mesh)
% trace_vertex traces the trajectory definied by Vertex vertex and velocity
% vector velo.
%
%  vertex is the foot of the trajectory, vertex_q and vertex_r the other two
%  vertices of this element, gradN is the gradient of basis function
%  associated with vertex
%
% Aloc(i,[1 ,2]) is head of trajectory
% Aloc(i,3) is the Element the head lies in
 

% u_normed is unit direction of the extrusion
u_length = norm(velo);
u_normed = velo / u_length;

% initial values for while loop

% tau is intersection parameter of extrusion and first edge
tau = -1/(gradN*u_normed');
% y is intersection point
y = vertex + tau * u_normed;
tau_total = tau;
% r and q are the other vertices in element i
r = vertex_r;
q = vertex_q;
Tid =Element;

% edge is global ID of the edge opposite to vertex
edge=Mesh.Vert2Edge(r_id,q_id);
% next triangle
h_Tid=Mesh.Edge2Elem(edge,:);
nTid=h_Tid(h_Tid~=Tid);

while (tau_total < u_length  && nTid~=0 )
    % the new vertex s
    vid_loc = Mesh.Elements(nTid,:);
    s_id=vid_loc(vid_loc~=r_id);
    s_id=s_id(s_id~=q_id);
    s = Mesh.Coordinates(s_id,:);

    % Compute element mapping
    BK_loc = [r-s; ...
        q-s];
    det_BK_loc = abs(det(BK_loc));
    inv_BK_loc = inv(BK_loc);

    % gradient of shape functions
    gradNr=[1,0]*inv_BK_loc';
    gradNq=[0,1]*inv_BK_loc';

    % barycentric coordinates of y
    Nr = norm(y-q)/norm(q-r); 
    Nq = norm(y-r)/norm(q-r); 

    % tau_q and tau_r are intersection parameters of extrusion and
    % edges opposite to r and q;
    
    % I) y==r and u_normed points into nTid
    f=-(gradNq*u_normed')/(gradNr*u_normed'); %
    % f is q-Coordinate of r+\tau u_normed  
    if (abs(Nr-1)<eps &&  f<=1 && f>=0)
        tau_r = -1/(gradNr*u_normed');
        y=y+tau_r*u_normed;
        tau_total=tau_total+tau_r;
        r_id=s_id;
        r=s;
    else
        %II) y==q and u_normed points into nTid
        f=-(gradNr*u_normed')/(gradNq*u_normed');
        % f is r-Coordinate of q+\tau u_normed
        if (abs(Nq-1)<eps && f-eps<=1 && f+eps>=0)
            tau_q = -1/(gradNq*u_normed');
            y=y+tau_q*u_normed;
            tau_total=tau_total+tau_q;
            q_id=s_id;
            q=s;
        else
            % III): not I) and not II)
            % we need to choose the minimal positive value if neither y==r
            % nor y==q;
            % In the case y==r (y==q) tau_q==0 (tau_r==0) and we
            % rotate around y until we end up in I) or II)
            tau_q = -Nq/(gradNq*u_normed');
            if (tau_q<0)  tau_q=inf; end
            tau_r = -Nr/(gradNr*u_normed');
            if (tau_r < 0)  tau_r = inf; end
            if (tau_r <=  tau_q)
                y=y+tau_r*u_normed;
                tau_total=tau_total+tau_r;
                r_id=s_id;
                r=s;
            end
            if (tau_q < tau_r)
                y=y+tau_q*u_normed;
                tau_total=tau_total+tau_q;
                q_id=s_id;
                q=s;
            end
        end
    end
    if (tau_q==inf && tau_r==inf)
        inf
    end
    % next triangle
    edge=Mesh.Vert2Edge(r_id,q_id);
    Tid=nTid;
    h_Tid=Mesh.Edge2Elem(edge,:);
    nTid=h_Tid(h_Tid~=Tid);
end %while
Aloc = [vertex+u_length*u_normed,Tid];