function varargout = assemMat_Lie_Inn_W1F(Mesh,V_HANDLE,varargin)

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang & Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  nElements = size(Mesh.Elements,1);
  nEdges = size(Mesh.Edges,1);
  
  % Preallocate memory
  
  I = zeros(2*9*nEdges,1);
  J = zeros(2*9*nEdges,1);
  A = zeros(2*9*nEdges,1);

  % Assemble element contributions
  
  loc = 1:9;
  for i = 1:nEdges
    
    Mloc = zeros(3,3);  
      
    % Extract vertices of current element
    
    vidx = Mesh.Edges(i,:);
    tangent = Mesh.Coordinates(vidx(2),:)-Mesh.Coordinates(vidx(1),:);
    mid=1/2*(Mesh.Coordinates(vidx(1),:)+Mesh.Coordinates(vidx(2),:));  
    Elem = Mesh.Edge2Elem(i,:);

    % Evaluate coefficient function at quadrature nodes
    Fval = 1/2*(V_HANDLE(Mesh.Coordinates(vidx(2),:),0,varargin{:})+...
        V_HANDLE(Mesh.Coordinates(vidx(1),:),0,varargin{:}));
    Fval =[Fval(:,2) -Fval(:,1)];
    
    signed_vol=(Fval*tangent');
    volume = abs(signed_vol);
     
    if abs(volume)>eps
    
        if (Elem(1)>0 && Elem(2)>0)
            % Right Element
            vidx_R=Mesh.Elements(Elem(2),:);
            Vertices_R = Mesh.Coordinates(vidx_R,:);

            P1_R = Vertices_R(1,:);
            P2_R = Vertices_R(2,:);
            P3_R = Vertices_R(3,:);

            % transpose of transformation matrix
            BK_R = [ P2_R - P1_R ; P3_R - P1_R ];
            det_BK_R = abs(det(BK_R));
            inv_BK_R = inv(BK_R);
            TK_R = transpose(inv_BK_R);

            % evaluation point
            mid_hat=(mid-P1_R)*inv_BK_R;
            N_R =shap_W1F(mid_hat);

            % Extract global edge numbers
            eidx_R = [Mesh.Vert2Edge(Mesh.Elements(Elem(2),2),Mesh.Elements(Elem(2),3)) ...
                Mesh.Vert2Edge(Mesh.Elements(Elem(2),3),Mesh.Elements(Elem(2),1)) ...
                Mesh.Vert2Edge(Mesh.Elements(Elem(2),1),Mesh.Elements(Elem(2),2))];

            % Determine the orientation

            if(Mesh.Edges(eidx_R(1),1)==vidx_R(2)),  p1_R = 1;  else    p1_R = -1;  end
            if(Mesh.Edges(eidx_R(2),1)==vidx_R(3)),  p2_R = 1;  else    p2_R = -1;  end
            if(Mesh.Edges(eidx_R(3),1)==vidx_R(1)),  p3_R = 1;  else    p3_R = -1;  end

            N_R([1 2]) = p1_R*N_R([1 2])*TK_R;
            N_R([3 4]) = p2_R*N_R([3 4])*TK_R;
            N_R([5 6]) = p3_R*N_R([5 6])*TK_R;

            % left Element
            vidx_L=Mesh.Elements(Elem(1),:);
            Vertices_L = Mesh.Coordinates(vidx_L,:);

            P1_L = Vertices_L(1,:);
            P2_L = Vertices_L(2,:);
            P3_L = Vertices_L(3,:);

            % transpose of transformation matrix
            BK_L = [ P2_L - P1_L ; P3_L - P1_L ];
            det_BK_L = abs(det(BK_L));
            inv_BK_L = inv(BK_L);
            TK_L = transpose(inv_BK_L);

            % evaluation point
            mid_hat=(mid-P1_L)*inv_BK_L;
            N_L =shap_W1F(mid_hat);

            % Extract global edge numbers
            eidx_L = [Mesh.Vert2Edge(Mesh.Elements(Elem(1),2),Mesh.Elements(Elem(1),3)) ...
                Mesh.Vert2Edge(Mesh.Elements(Elem(1),3),Mesh.Elements(Elem(1),1)) ...
                Mesh.Vert2Edge(Mesh.Elements(Elem(1),1),Mesh.Elements(Elem(1),2))];

            % Determine the orientation
            if(Mesh.Edges(eidx_L(1),1)==vidx_L(2)),  p1_L = 1;  else    p1_L = -1;  end
            if(Mesh.Edges(eidx_L(2),1)==vidx_L(3)),  p2_L = 1;  else    p2_L = -1;  end
            if(Mesh.Edges(eidx_L(3),1)==vidx_L(1)),  p3_L = 1;  else    p3_L = -1;  end

            N_L([1 2]) = p1_L*N_L([1 2])*TK_L;
            N_L([3 4]) = p2_L*N_L([3 4])*TK_L;
            N_L([5 6]) = p3_L*N_L([5 6])*TK_L;

            % Compute local mass matrix
            Mloc(1,1) = sum(volume.*sum(N_L([1 2]).*N_R([1 2]),2));
            Mloc(2,1) = sum(volume.*sum(N_L([3 4]).*N_R([1 2]),2));
            Mloc(3,1) = sum(volume.*sum(N_L([5 6]).*N_R([1 2]),2));

            Mloc(1,2) = sum(volume.*sum(N_L([1 2]).*N_R([3 4]),2));
            Mloc(2,2) = sum(volume.*sum(N_L([3 4]).*N_R([3 4]),2));
            Mloc(3,2) = sum(volume.*sum(N_L([5 6]).*N_R([3 4]),2));

            Mloc(1,3) = sum(volume.*sum(N_L([1 2]).*N_R([5 6]),2));
            Mloc(2,3) = sum(volume.*sum(N_L([3 4]).*N_R([5 6]),2));
            Mloc(3,3) = sum(volume.*sum(N_L([5 6]).*N_R([5 6]),2));

            % Add contributions to stiffness matrix

            if signed_vol>0 % left is downwind

                % match orientation to left element
                Eloc = Mesh.EdgeLoc(i,2);
                if(vidx(1) == Mesh.Elements(Elem(2),rem(Eloc,3)+1))
                    Match = 1;
                else
                    Match = -1;
                end

                Mloc=Mloc;

                I(loc) = set_Rows(eidx_L,3);
                J(loc) = set_Cols(eidx_R,3);
                A(loc) = - Mloc(:);
                loc = loc+9;

                % Compute local mass matrix

                Mloc_L(1,1) = sum(volume.*sum(N_L([1 2]).*N_L([1 2]),2));
                Mloc_L(2,1) = sum(volume.*sum(N_L([3 4]).*N_L([1 2]),2));
                Mloc_L(3,1) = sum(volume.*sum(N_L([5 6]).*N_L([1 2]),2));

                Mloc_L(1,2) = sum(volume.*sum(N_L([1 2]).*N_L([3 4]),2));
                Mloc_L(2,2) = sum(volume.*sum(N_L([3 4]).*N_L([3 4]),2));
                Mloc_L(3,2) = sum(volume.*sum(N_L([5 6]).*N_L([3 4]),2));

                Mloc_L(1,3) = sum(volume.*sum(N_L([1 2]).*N_L([5 6]),2));
                Mloc_L(2,3) = sum(volume.*sum(N_L([3 4]).*N_L([5 6]),2));
                Mloc_L(3,3) = sum(volume.*sum(N_L([5 6]).*N_L([5 6]),2));

                I(loc) = set_Rows(eidx_L,3);
                J(loc) = set_Cols(eidx_L,3);
                A(loc) = Mloc_L(:);
                loc = loc+9;

            end

            if signed_vol<0 % right is downwind

                % match orientation to right element
                Eloc = Mesh.EdgeLoc(i,1);
                if(vidx(1) == Mesh.Elements(Elem(1),rem(Eloc,3)+1))
                    Match = 1;
                else
                    Match = -1;
                end

                Mloc=Mloc';

                I(loc) = set_Rows(eidx_R,3);
                J(loc) = set_Cols(eidx_L,3);
                A(loc) = -Mloc(:);
                loc = loc+9;

                % Compute local mass matrix

                Mloc_R(1,1) = sum(volume.*sum(N_R([1 2]).*N_R([1 2]),2));
                Mloc_R(2,1) = sum(volume.*sum(N_R([3 4]).*N_R([1 2]),2));
                Mloc_R(3,1) = sum(volume.*sum(N_R([5 6]).*N_R([1 2]),2));

                Mloc_R(1,2) = sum(volume.*sum(N_R([1 2]).*N_R([3 4]),2));
                Mloc_R(2,2) = sum(volume.*sum(N_R([3 4]).*N_R([3 4]),2));
                Mloc_R(3,2) = sum(volume.*sum(N_R([5 6]).*N_R([3 4]),2));

                Mloc_R(1,3) = sum(volume.*sum(N_R([1 2]).*N_R([5 6]),2));
                Mloc_R(2,3) = sum(volume.*sum(N_R([3 4]).*N_R([5 6]),2));
                Mloc_R(3,3) = sum(volume.*sum(N_R([5 6]).*N_R([5 6]),2));

                I(loc) = set_Rows(eidx_R,3);
                J(loc) = set_Cols(eidx_R,3);
                A(loc) = Mloc_R(:);
                loc = loc+9;

            end

        end
    end
  end
  % Assign output arguments
  
  I=I(1:loc(1)-1);
  J=J(1:loc(1)-1);
  A=A(1:loc(1)-1);
  
  if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
  else
    varargout{1} = sparse(I,J,A,nEdges,nEdges);      
  end
  
return