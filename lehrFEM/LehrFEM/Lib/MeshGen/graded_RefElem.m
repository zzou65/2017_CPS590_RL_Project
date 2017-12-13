function Mesh = graded_RefElem(nNodes,FHandle,varargin)
% GRADED_REFELEM Generates a graded point set inside the reference element. 
%
%   MESH = GRADED_REFELEM(NNODES,FHANDLE) generates the points set of the
%   graded mesh in the reference element by the given number of nodes and
%   the algebraically graded function.
%
%   Example:
%
%   Mesh = graded_RefElem(10,FHandle);

%   Copyright 2005-2005 Patrick Meury & Mengyu Wang
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

    % Preallocate memory
    
    Label = zeros(nNodes+1,1);
    DEdge_Length = zeros(nNodes,1);
    AEdge_Length = zeros(nNodes,1);
    nTNodes = zeros(nNodes,1);
    nCNoees = zeros(nNodes,1);
    
    % Extract information from algebraically graded division
    
    Label = FHandle(linspace(0,1,nNodes+1),varargin{:});
    Width = Label(2:nNodes+1) - Label(1:nNodes);
    AEdge_Length = 2^.5*Label(2:nNodes+1);
    
    % Compute the number of nodes in each layer
    
    nTNodes = genNodes(nNodes,Width,AEdge_Length);
    nCNodes = cumsum(nTNodes)+1;
    
    % Generate Coordinates field
    
    Mesh.Coordinates = zeros(1+sum(nTNodes),2);

    for i = 1:nNodes
        Para = linspace(0,1,nTNodes(i));
        if(i==1)
            Mesh.Coordinates(2:nCNodes(i),:) = Label(i+1)*[(1-Para') Para'];
        else
            Mesh.Coordinates(1+nCNodes(i-1):nCNodes(i),:) = Label(i+1)*[(1-Para') Para'];
        end
    end

return

%%%%%%%%%%%%% Node generating strategy %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function nTNodes = genNodes(nNodes,width,length)

  % Preallocate memory
  
  nTNodes = zeros(nNodes,1);
  
  % Compute nTNodes
  
  nTNodes = (1:nNodes)'+1;

return
