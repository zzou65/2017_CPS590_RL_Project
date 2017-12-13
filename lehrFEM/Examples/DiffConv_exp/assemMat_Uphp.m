function varargout = assemMat_Uphp(Mesh, Elem2Dof, VHandle,QuadRule, Shap,pmax)
% assemMat_UpHP  
%    Up: Upwind part, using DOFS on element boundary
%
%   Copyright 2008-2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

 % Initialize constants
 % Assign output arguments
 nF=size(Mesh.Elements,1);
 nE=size(Mesh.Edges,1);
 nC=size(Mesh.Coordinates,1);
  
 % Preallocate memory
  
 % edge and vertices sourounding volume
 ME = zeros(nE,1);
 MV = zeros(nC,1);
  
 for f = 1:nF
     % Vertices
     vid = Mesh.Elements(f,:);
     a1 = Mesh.Coordinates(vid(1),:);
     a2 = Mesh.Coordinates(vid(2),:);
     a3 = Mesh.Coordinates(vid(3),:);
     
     % Element Mass
     BK = [a2-a1;a3-a1];
     det_BK = abs(det(BK));
     
     % Extract global edge numbers
     eidx =... 
           [Mesh.Vert2Edge(Mesh.Elements(f,1),Mesh.Elements(f,2)) ...
            Mesh.Vert2Edge(Mesh.Elements(f,2),Mesh.Elements(f,3)) ...
            Mesh.Vert2Edge(Mesh.Elements(f,3),Mesh.Elements(f,1))];
     
     ME(eidx) = ME(eidx) + det_BK*[1 ;1; 1];
     MV(vid) = MV(vid)+ det_BK*[1 ;1; 1];
 end

 % upwinding part
 
 % Preallocate memory
    
  nEntries = sum((3*ones(nF,1) + Elem2Dof.EDofs{1}.nDofs + ...
                  Elem2Dof.EDofs{2}.nDofs + Elem2Dof.EDofs{3}.nDofs + ...
                  Elem2Dof.CDofs.nDofs).^2);
  I = zeros(nEntries,1);
  J = zeros(nEntries,1);
  A = zeros(nEntries,1);
    
  % Assemble global stiffness matrix  
    
  EDofs = zeros(1,3);
  EDir = zeros(1,3);
  CDofs = 0;
  offset = 0;
 
  for i = 1:nF
  
    % Extract vertices of current element
    
    vidx = Mesh.Elements(i,:);
    Vertices = Mesh.Coordinates(vidx,:);
    eidx = [...
           Mesh.Vert2Edge(Mesh.Elements(i,2),Mesh.Elements(i,3))...
           Mesh.Vert2Edge(Mesh.Elements(i,3),Mesh.Elements(i,1))...
           Mesh.Vert2Edge(Mesh.Elements(i,1),Mesh.Elements(i,2))];  
    
    % Extract local polynomial orders
    
    EDofs(1) = Elem2Dof.EDofs{1}.nDofs(i);
    EDofs(2) = Elem2Dof.EDofs{2}.nDofs(i);
    EDofs(3) = Elem2Dof.EDofs{3}.nDofs(i); 
    CDofs = Elem2Dof.CDofs.nDofs(i);
      
    % Extract local edge orientations
    
    EDir(1) = Elem2Dof.EDofs{1}.Dir(i);
    EDir(2) = Elem2Dof.EDofs{2}.Dir(i);
    EDir(3) = Elem2Dof.EDofs{3}.Dir(i);
    
    % Extract local volume contributions
    
    EMass = ME(eidx); 
    VMass = MV(vidx);
    
    % Compute element contributions
    
    Aloc = STIMA_UPhp(Vertices, Mesh.ElemFlag(i), EDofs, EDir,...
        CDofs, EMass, VMass, QuadRule, VHandle, Shap,pmax);
    
    % Add contributions to global matrix
    
    idx = [vidx ...
           Elem2Dof.EDofs{1}.Dofs{i} ...
           Elem2Dof.EDofs{2}.Dofs{i} ...
           Elem2Dof.EDofs{3}.Dofs{i} ...
           Elem2Dof.CDofs.Dofs{i}];
    n_idx = 3+sum(EDofs)+CDofs;
    
    loc = offset + (1:n_idx^2); 
    I(loc) = set_Rows(idx,n_idx);
    J(loc) = set_Cols(idx,n_idx);
    A(loc) = Aloc(:);
    
    offset = offset + n_idx^2;
    
  end  
      
  % Assign output arguments
  
  if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = A;
  else
    varargout{1} = sparse(I,J,A);    
  end;