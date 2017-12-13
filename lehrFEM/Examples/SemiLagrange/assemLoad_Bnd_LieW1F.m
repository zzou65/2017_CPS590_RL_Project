function L = assemLoad_Bnd_LieW1F(Mesh,BdFlags,L,QuadRule,FHandle, VHandle,varargin)
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    EDGES        P-by-2 matrix specifying the edges of the mesh.
%    BDFLAGS      P-by-1 matrix specifying the boundary type of each boundary
%                 edge in the mesh.
%    EDGE2ELEM    N-by-2 matrix connecting edges to elements. The first column
%                 specifies the left hand side element where the second column
%                 specifies the right hand side element.
%    EDGELOC      P-by-2 matrix connecting egdes to local edges of elements.
%
%   See also GET_BDEDGES, SHAP_LFE.

%   Copyright 2005-2005 Patrick Meury & Kah-Ling Sia
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constants

nCoordinates = size(Mesh.Coordinates,1);
nGauss = size(QuadRule.w,1);

% Precompute shape functions

Lloc = zeros(3,1);
for j1 = BdFlags

    % Extract Neumann edges

    Loc = get_BdEdges(Mesh);
    Loc = Loc(Mesh.BdFlags(Loc) == j1);

    for j2 = Loc'

        % Compute element map
        
        Edge = Mesh.Coordinates(Mesh.Edges(j2,:),:);
        Normal = Mesh.Normals(j2,:);

        % Extract left or right hand side element data

        if(Mesh.Edge2Elem(j2,1) > 0)
            Data.Element = Mesh.Edge2Elem(j2,1);
%            Data.ElemFlag = ElemFlag(Data.Element);
            Data.Vertices = Mesh.Coordinates(Mesh.Elements(Data.Element,:),:);
            Data.EdgeLoc = Mesh.EdgeLoc(j2,1);
            Data.Match = Mesh.EdgeOrient(j2,1);
        else
            Data.Element = Mesh.Edge2Elem(j2,2);
%            Data.ElemFlag = ElemFlag(Data.Element);
            Data.Vertices = Mesh.Coordinates(Mesh.Elements(Data.Element,:),:);
            Data.EdgeLoc = Mesh.EdgeLoc(j2,2);
            Data.Match = Mesh.EdgeOrient(j2,2);
        end
        
        % Compute element contributions
        nPoints = size(QuadRule.w,1);

        dS = norm(Edge(2,:)-Edge(1,:));
        x = QuadRule.x*(Edge(2,:)-Edge(1,:)) + ones(nPoints,1)*Edge(1,:);
        Vval= VHandle(x);
        Fval= FHandle(x,varargin{:});

        % Compute on the element
        
        bK = Data.Vertices(1,:);
        BK = [Data.Vertices(2,:)-bK; ...
            Data.Vertices(3,:)-bK];
        inv_BK_t = transpose(inv(BK));
        det_BK=abs(det(BK));

        % Compute constant gradients of barycentric coordinate functions
        g1 = [Data.Vertices(2,2)-Data.Vertices(3,2);Data.Vertices(3,1)-Data.Vertices(2,1)]/det_BK;
        g2 = [Data.Vertices(3,2)-Data.Vertices(1,2);Data.Vertices(1,1)-Data.Vertices(3,1)]/det_BK;
        g3 = [Data.Vertices(1,2)-Data.Vertices(2,2);Data.Vertices(2,1)-Data.Vertices(1,1)]/det_BK;

        % Get barycentric coordinates of quadrature points
        if Data.Match
            switch(Data.EdgeLoc)
                case 1
                    baryc = [zeros(nPoints,1) 1-QuadRule.x QuadRule.x];
                case 2
                    baryc = [QuadRule.x zeros(nPoints,1) 1-QuadRule.x];
                case 3
                    baryc = [1-QuadRule.x QuadRule.x zeros(nPoints,1)];
            end
        else
            switch(Data.EdgeLoc)
                case 1
                    baryc = [zeros(nPoints,1) QuadRule.x 1-QuadRule.x];
                case 2
                    baryc = [1-QuadRule.x zeros(nPoints,1) QuadRule.x];
                case 3
                    baryc = [QuadRule.x 1-QuadRule.x zeros(nPoints,1)];
            end
        end
        N1 = baryc(:,2)*g3'-baryc(:,3)*g2';
        N2= baryc(:,3)*g1'-baryc(:,1)*g3';
        N3 = baryc(:,1)*g2'-baryc(:,2)*g1';

        N_F=Fval*Normal';
        
        N_V=Vval*Normal';
        
        V_F=sum(Vval.*Fval,2);
        
        N_N1=N1*Normal';
        N_N2=N2*Normal';
        N_N3=N3*Normal';

        % Compute entries of element matrix

        %   if(LData.Element < RData.Element)
        %     gamma = 1;
        %   else
        %     gamma = -1;
        %   end
        gamma=dS;
%         Lloc(1) = gamma*sum(QuadRule.w.*N_N1.*N_F.*N_V);
%         Lloc(2) = gamma*sum(QuadRule.w.*N_N2.*N_F.*N_V);
%         Lloc(3) = gamma*sum(QuadRule.w.*N_N3.*N_F.*N_V);
        
        Lloc(1) = gamma*sum(QuadRule.w.*N_N1.*V_F);
        Lloc(2) = gamma*sum(QuadRule.w.*N_N2.*V_F);
        Lloc(3) = gamma*sum(QuadRule.w.*N_N3.*V_F);

        % Add element contributions to stiffness matrix
        idx = [Mesh.Vert2Edge(Mesh.Elements(Data.Element,2),Mesh.Elements(Data.Element,3)) ...
            Mesh.Vert2Edge(Mesh.Elements(Data.Element,3),Mesh.Elements(Data.Element,1)) ...
            Mesh.Vert2Edge(Mesh.Elements(Data.Element,1),Mesh.Elements(Data.Element,2))];

        % Determine the orientation
        if(Mesh.Edges(idx(1),1)==Mesh.Elements(Data.Element,2)),  p1 = 1;  else    p1 = -1;  end
        if(Mesh.Edges(idx(2),1)==Mesh.Elements(Data.Element,3)),  p2 = 1;  else    p2 = -1;  end
        if(Mesh.Edges(idx(3),1)==Mesh.Elements(Data.Element,1)),  p3 = 1;  else    p3 = -1;  end

        Peori = [p1; p2; p3]; % scaling matrix taking into account orientations
        Lloc = Peori.*Lloc;

        L(idx) = L(idx)+Lloc;
    end
end

return