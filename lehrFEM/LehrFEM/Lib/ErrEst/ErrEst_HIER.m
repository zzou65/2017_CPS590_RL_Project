function Eta = ErrEst_HIER(U,Mesh,RHandle,SHandle,varargin)
% ERREST_HIER Hierarchical error estimator.
%
%   ETA = ERREST_HIER(U,MESH,RHANDLE,SHANDLE) computes the hierarchical 
%   error estimator ETA of the solution U.
%   
%   The struct MESH must at least contain the following fields:
%    COORDINATES  M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS     N-by-3 matrix specifying the elements of the mesh.
%    EDGE2ELEM    P-by-2 matrix connecting edges to elements. The first 
%                 column specifies the left hand side element where the
%                 second column specifies the right hand side element.
%    EDGELOC      P-by-3 matrix connecting egdes to local edges of
%                 elements.
%
%   RHANDLE is function handle to the element weak residual and SHANDLE is
%   function handle to the diagonal element stiffness matrix. 
%
%   ETA = ERREST_HIER(U,MESH,RHANDLE,SHANDLE,PARAM) also handles the
%   variable length argument list PARAM to the element residual RHANDLE and
%   the element stiffness matrix SHANDLE,
%
%   Example:
%
%   Eta = ErrEst_HIER(U,Mesh,@Res_Lapl,@STIMA_Lapl_ErrEst);

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland
  
  % Initialize constants
   
  nElements = size(Mesh.Elements,1);
  nEdges = size(Mesh.Edges,1);
  
  % Preallocate memory
  
  Eta = zeros(nElements,1);        
  for i = 1:nEdges
    if(Mesh.BdFlags(i) >= 0)
    
      % Extract vertices of left and right hand side elements
      
      Elem_l = Mesh.Edge2Elem(i,1);
      Elem_r = Mesh.Edge2Elem(i,2);      
      vidx_l = Mesh.Elements(Elem_l,:);
      vidx_r = Mesh.Elements(Elem_r,:);
      
      % Compute left and right hand-side residual for the current edge
      
      Res = RHandle(U(vidx_l), ... 
                    Mesh.EdgeLoc(i,1), ...
                    Mesh.Coordinates(vidx_l,:), ...
                    Mesh.ElemFlag(Elem_l), ...
                    varargin{:}) + ...
            RHandle(U(vidx_r), ...
                    Mesh.EdgeLoc(i,2), ...
                    Mesh.Coordinates(vidx_r,:), ...
                    Mesh.ElemFlag(Elem_r), ...
                    varargin{:});
      
      % Compute diagonal entry of the global stiffness matrix on the current edge          
                
      Diag_l = SHandle(Mesh.EdgeLoc(i,1), ...
                       Mesh.Coordinates(vidx_l,:), ...
                       Mesh.ElemFlag(Elem_l), ...
                       varargin{:});
      Diag_r = SHandle(Mesh.EdgeLoc(i,2), ...
                       Mesh.Coordinates(vidx_r,:), ...
                       Mesh.ElemFlag(Elem_r), ...
                       varargin{:});
                 
      % Add edge contributions to left and right hand side neighbours
      
      ECoeff = Res/(Diag_l+Diag_r);
      Eta(Elem_l) = Eta(Elem_l) + abs(ECoeff*Diag_l*ECoeff);
      Eta(Elem_r) = Eta(Elem_r) + abs(ECoeff*Diag_r*ECoeff);
                 
    end
  end
  Eta = sqrt(Eta);
  
return