function varargout = assemMat_SemiLag_W1Fcheat(Mesh, tracedvertices, varargin)
% assemMat_SemiLag_W1F calculates the interpolation matrix for
% Semi-Lagrange for Whitney 1-forms.
%
% Mesh structure; tracedvertices array with coodinates of traced vertices and
% local local elements.
%
% Example:
%
%   Copyright 2008-2009 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

nEdges =size(Mesh.Edges,1);

% Preallocate memory
pbE = [];
I=[];
J=[];
level=[1 2 3];

% loop over all edges
for i = 1:nEdges
    %i
    %  if i==1811
    %        i
    % end
    % endpoints
    epid = Mesh.Edges(i,:);
    e1 = tracedvertices(epid(1),[1 2]);
    e2 = tracedvertices(epid(2),[1 2]);

    % Element for e1 and element mapping
    El_ID = tracedvertices(epid(1),3);
    
    if El_ID==0 && tracedvertices(epid(2),3)==0
        pbE(level(1)) = 1/3;
        pbE(level(2)) = 1/3;
        pbE(level(3)) = 1/3;
        I(level)=[i,i,i];
        J(level)=[i,i,i];
        level=level+3;
    else
    	if El_ID == 0
 	       e1 = tracedvertices(epid(2),[1 2]);
               e2 = tracedvertices(epid(1),[1 2]);
	       El_ID = tracedvertices(epid(2),3);
        end

        vid = Mesh.Elements(El_ID,:);
        a1 = Mesh.Coordinates(vid(1),:);
        a2 = Mesh.Coordinates(vid(2),:);
        a3 = Mesh.Coordinates(vid(3),:);

        % Extract global edge numbers
        eidx = [Mesh.Vert2Edge(vid(2),vid(3)) ...
            Mesh.Vert2Edge(vid(3),vid(1)) ...
            Mesh.Vert2Edge(vid(1),vid(2))];

        bK = a1;
        BK = [a2-bK; ...
            a3-bK];
        inv_BK = inv(BK);

        % gradient of shape functions
        gradN1 = [-1,-1]*inv_BK';
        gradN2 = [1,0]*inv_BK';
        gradN3 = [0,1]*inv_BK';

        % edge direction
        u = e2-e1;
        u_length=norm(u);
        u_normed=u/u_length;

        % initial values for while loop
        tau=0;
        edge=0;

        % tau is intersection parameter in direction u_normed with edges
        e1_hat=(e1-bK)*inv_BK;

        if (norm(e1-a1)<eps || norm(e1-a2)<eps || norm(e1-a3)<eps)
            % I) e1 is a vertex
            % find Element u_normed points into

            El_IDs=nonzeros(Mesh.AdjElements(epid(1),:));
            counter=1;
            notfound=1;

            while (notfound)
                El_ID=El_IDs(counter);

                vid = Mesh.Elements(El_ID,:);
                a1 = Mesh.Coordinates(vid(1),:);
                a2 = Mesh.Coordinates(vid(2),:);
                a3 = Mesh.Coordinates(vid(3),:);

                %Extract global edge numbers
                eidx = [Mesh.Vert2Edge(vid(2),vid(3)) ...
                    Mesh.Vert2Edge(vid(3),vid(1)) ...
                    Mesh.Vert2Edge(vid(1),vid(2))];

                bK = a1;
                BK = [a2-bK; ...
                    a3-bK];
                inv_BK = inv(BK);

                %gradient of shape functions
                gradN1 = [-1,-1]*inv_BK';
                gradN2 = [1,0]*inv_BK';
                gradN3 = [0,1]*inv_BK';

                if (epid(1)==vid(1) && abs(gradN1*u_normed')>eps)
                    tau = -1/(gradN1*u_normed');
                    f = tau*(gradN2*u_normed');
                    edge = eidx(1);
                end
                if (epid(1)==vid(2) && abs(gradN2*u_normed')>eps)
                    tau = -1/(gradN2*u_normed');
                    f = tau*(gradN3*u_normed');
                    edge = eidx(2);
                end
                if (epid(1) == vid(3) && abs(gradN3*u_normed')>eps)
                    tau = -1/(gradN3*u_normed');
                    f = tau*(gradN1*u_normed');
                    edge = eidx(3);
                end
                % HH if (f-eps<=1 && f+eps>=0 && tau>0) notfound=0; end
                if (f-1<=eps && f+100*eps>=0 && tau>0) notfound=0; end
                counter = counter + 1;
            end

        else
            if (abs(sum(e1_hat)-1)<=eps || abs(e1_hat(1))<=eps || abs(e1_hat(2))<=eps)
                % II) e1 is on a edge
                % find element u_normed points into
                if abs(sum(e1_hat)-1)<=eps
                    or = (a3-a2)*[0 1; -1 0]*u_normed';
                    c_edge=eidx(1);
                else if abs(e1_hat(1))<=eps
                        or = (a1-a3)*[0 1; -1 0]*u_normed';
                        c_edge=eidx(2);
                    else
                        or = (a2-a1)*[0 1; -1 0]*u_normed';
                        c_edge=eidx(3);
                    end
                end

                if (sign(or) == -1)
                    h_Tid=Mesh.Edge2Elem(c_edge,:);
                    El_ID=h_Tid(h_Tid~=El_ID);

                    vid = Mesh.Elements(El_ID,:);
                    a1 = Mesh.Coordinates(vid(1),:);
                    a2 = Mesh.Coordinates(vid(2),:);
                    a3 = Mesh.Coordinates(vid(3),:);

                    % Extract global edge numbers
                    eidx = [Mesh.Vert2Edge(vid(2),vid(3)) ...
                        Mesh.Vert2Edge(vid(3),vid(1)) ...
                        Mesh.Vert2Edge(vid(1),vid(2))];

                    bK = a1;
                    BK = [a2-bK; ...
                        a3-bK];
                    inv_BK = inv(BK);

                    % gradient of shape functions
                    gradN1 = [-1,-1]*inv_BK';
                    gradN2 = [1,0]*inv_BK';
                    gradN3 = [0,1]*inv_BK';

                    e1_hat=(e1-bK)*inv_BK;
                end
            end
            % e1 is inside the element
            tau=zeros(3,1);
            h=gradN1*u_normed';
            if (abs(h)>eps)  tau(1) = -(1-sum(e1_hat))/h;  end
            h=gradN2*u_normed';
            if (abs(h)>eps)  tau(2) = -e1_hat(1)/h;  end
            h=gradN3*u_normed';
            if (abs(h)>eps)  tau(3) =  -e1_hat(2)/h;  end

            tau(tau<=0)=inf;
            [tau, edge_loc]=min(tau);
            edge=eidx(edge_loc);
        end

        % next triangle Tid
        h_Tid=Mesh.Edge2Elem(edge,:);
        Tid=h_Tid(h_Tid~=El_ID);

        % y is intersection point
        x = e1;
        y = e1 + tau * u_normed;
        tau_total = tau;
        % r and q are the other vertices of  edge i
        r_id = Mesh.Edges(edge,1);
        q_id = Mesh.Edges(edge,2);
        r = Mesh.Coordinates(r_id,:);
        q = Mesh.Coordinates(q_id,:);
        calculate=1;

        while (tau_total < u_length  && Tid~=0 )

            if calculate
                %  contribution of Element El_ID

                vid = Mesh.Elements(El_ID,:);
                a1 = Mesh.Coordinates(vid(1),:);
                a2 = Mesh.Coordinates(vid(2),:);
                a3 = Mesh.Coordinates(vid(3),:);

                % Extract global edge numbers
                eidx = [Mesh.Vert2Edge(vid(2),vid(3)) ...
                    Mesh.Vert2Edge(vid(3),vid(1)) ...
                    Mesh.Vert2Edge(vid(1),vid(2))];

                bK = a1;
                BK = [a2-bK; ...
                    a3-bK];
                inv_BK = inv(BK);

                % compute contribution of edges of Element El_ID
                x_hat=(x-bK)*inv_BK;
                y_hat=(y-bK)*inv_BK;

                % Determine the orientation
                if(Mesh.Edges(eidx(1),1)==vid(2)),  p1 = 1;  else    p1 = -1;  end
                if(Mesh.Edges(eidx(2),1)==vid(3)),  p2 = 1;  else    p2 = -1;  end
                if(Mesh.Edges(eidx(3),1)==vid(1)),  p3 = 1;  else    p3 = -1;  end

                % Pathintegrals for edge basis functions
                contrib = x_hat*[0 1; -1 0]*y_hat';
                pbE(level(1)) = p1 * (contrib);
                pbE(level(2)) = p2 * (contrib+x_hat(2)-y_hat(2));
                pbE(level(3)) = p3 * (contrib-x_hat(1)+y_hat(1));

                I(level)=[i,i,i];
                J(level)=eidx;
                level=level+3;

            end

            % the new vertex s
            vid_loc = Mesh.Elements(Tid,:);
            s_id=vid_loc(vid_loc~=r_id);
            s_id=s_id(s_id~=q_id);
            s = Mesh.Coordinates(s_id,:);

            % Compute element mapping
            bK=s;
            BK = [r-bK; ...
                q-bK];
            inv_BK = inv(BK);

            % gradient of shape functionsq
            gradNr=[1,0]*inv_BK';
            gradNq=[0,1]*inv_BK';

            % tau_q and tau_r are intersection parameters of extrusion and
            % edges opposite to r and q;
            % we need to choose the minimal positive value

            Nr = norm(y-q)/norm(q-r); % r-coordinate of y
            Nq = norm(y-r)/norm(q-r); % q-coordinate of y
            calculate = 0;
            if (norm(gradNr*u_normed')>eps)
                f=-(gradNq*u_normed')/(gradNr*u_normed');
                % f is q-Coordinate of r+\tau u_normed
            else
                f=inf;
            end
            % I) y==r and u_normed points into nTid
            if (abs(Nr-1)<eps &&  f-eps<=1 && f+eps>=0)
                tau_r = -1/(gradNr*u_normed');
                x=y;
                y=y+tau_r*u_normed;
                tau_total=tau_total+tau_r;
                r_id=s_id;
                r=s;
                calculate = 1;
            else
                if (norm(gradNq*u_normed') >eps )
                    f=-(gradNr*u_normed')/(gradNq*u_normed');
                else
                    f=inf;
                end
                % II) y==q and u_normed points into nTid
                if (abs(Nq-1)<eps && f-eps<=1 && f+eps>=0) %y is at vertex q and
                    tau_q = -1/(gradNq*u_normed');
                    x=y;
                    y=y+tau_q*u_normed;
                    tau_total=tau_total+tau_q;
                    q_id=s_id;
                    q=s;
                    calculate = 1;
                else
                    % III) not I) and not II)
                    % we need to choose the minimal positive value if neither y==r
                    % nor y==q;
                    % In the case y==r (y==q) tau_q==0 (tau_r==0) and we
                    % rotate around y until we end up in I) or II)
                    tau_q = inf;
                    if ((gradNq*u_normed')<0)  tau_q = -Nq/(gradNq*u_normed'); end
                    tau_r = inf;
                    if (gradNr*u_normed' < 0)  tau_r = -Nr/(gradNr*u_normed'); end
                    if (tau_r <=  tau_q)
                        if (abs(Nr-1)>eps) calculate = 1; end
                        x=y;
                        y = y  + tau_r*u_normed;
                        tau_total = tau_total + tau_r;
                        r_id = s_id;
                        r = s;
                    end
                    if (tau_q < tau_r)
                        if (abs(Nr-1)>eps) calculate = 1; end
                        x=y;
                        y=y+tau_q*u_normed;
                        tau_total=tau_total+tau_q;
                        q_id=s_id;
                        q=s;
                    end
                end
            end
            if (tau_q==inf && tau_r==inf)
                [i  inf];
            end

            % next triangle
            edge=Mesh.Vert2Edge(r_id,q_id);
            El_ID=Tid;
            h_Tid=Mesh.Edge2Elem(edge,:);
            Tid=h_Tid(h_Tid~=Tid);
            %Tid = setdiff(Mesh.Edge2Elem(edge,:),Tid);

        end %while
        % compute contribution of edges of Element El_ID
        y = e1 + u_length * u_normed;

        vid = Mesh.Elements(El_ID,:);
        a1 = Mesh.Coordinates(vid(1),:);
        a2 = Mesh.Coordinates(vid(2),:);
        a3 = Mesh.Coordinates(vid(3),:);

        % Extract global edge numbers
        eidx = [Mesh.Vert2Edge(vid(2),vid(3)) ...
            Mesh.Vert2Edge(vid(3),vid(1)) ...
            Mesh.Vert2Edge(vid(1),vid(2))];

        bK = a1;
        BK = [a2-bK; ...
            a3-bK];
        inv_BK = inv(BK);

        x_hat=(x-bK)*inv_BK;
        y_hat=(y-bK)*inv_BK;

        % Determine the orientation

        if(Mesh.Edges(eidx(1),1)==vid(2)),  p1 = 1;  else    p1 = -1;  end
        if(Mesh.Edges(eidx(2),1)==vid(3)),  p2 = 1;  else    p2 = -1;  end
        if(Mesh.Edges(eidx(3),1)==vid(1)),  p3 = 1;  else    p3 = -1;  end

        % Pathintegrals for edge basis functions
        contrib = x_hat*[0 1; -1 0]*y_hat';
        pbE(level(1)) = p1 * (contrib);
        pbE(level(2)) = p2 * (contrib+x_hat(2)-y_hat(2));
        pbE(level(3)) = p3 * (contrib-x_hat(1)+y_hat(1));
        I(level)=[i,i,i];
        J(level)=eidx;
        level=level+3;
    
    end %El_ID==0

end% for
if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = pbE;
else
    varargout{1} = sparse(I,J,pbE,nEdges,nEdges);
end
