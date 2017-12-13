function L = assemCochainD_1f(Mesh, FHandle, QuadRule_1D)
% assemCochain_1f assemble Cochain of one forms on dual grid.
%
%   L = ASSEMCochain_1f(MESH,FHandel) .... 
%
%   Example:
%
%   Mesh = load_Mesh('Coord_LShap.dat','Elem_LShap.dat');
%   A = assemMat_MassTwoD(Mesh);
%  
%   Copyright 2007-2007 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constants
% Assign output arguments
 
 nEdges=size(Mesh.Edges,1);
 L = zeros(nEdges,1);
 nGuass = size(QuadRule_1D.w,1);
 
 for i = 1:nEdges

        P1 = Mesh.Coordinates(Mesh.Edges(i,1),:);
        P2 = Mesh.Coordinates(Mesh.Edges(i,2),:);
        m=(P1+P2)/2;
        Elem=Mesh.Edge2Elem(i,:);
        if (Elem(1)~=0)
            vid=Mesh.Elements(Elem(1),:);
            P3=setdiff(vid,Mesh.Edges(i,:));
            b=(P1+P2+P3)/3;
            x = ones(nGuass,1)*m+QuadRule_1D.x*(b-m);
            dS = ones(nGuass,1)*(b-m);
            dS=[dS(:,2) dS(:,1)];
            Fval = FHandle(x);
            if (Elem(2)==0)
             L(i) = L(i)+sum(QuadRule_1D.w.*sum(Fval.*dS,2));
            else L(i) = L(i)+1/2sum(QuadRule_1D.w.*sum(Fval.*dS,2));
            end
        end
        if (Elem(2)~=0)
            vid=Mesh.Elements(Elem(2),:);
            P3=setdiff(vid,Mesh.Edges(i,:));
            b=(P1+P2+P3)/3;
            x = ones(nGuass,1)*m+QuadRule_1D.x*(b-m);
            dS = ones(nGuass,1)*(b-m);
            dS=[dS(:,2) dS(:,1)];
            Fval = FHandle(x);
            if (Elem(1)==0)
             L(i) = L(i)+sum(QuadRule_1D.w.*sum(Fval.*dS,2));
            else L(i) = L(i)+1/2sum(QuadRule_1D.w.*sum(Fval.*dS,2));
            end
        end
 end
  
return
   
  