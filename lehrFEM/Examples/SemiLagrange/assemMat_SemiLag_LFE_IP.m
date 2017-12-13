function varargout = assemMat_SemiLag_LFE_IP(Mesh, tracedvertices, varargin)
% assemMat_SemiLag_LFE_IP calculates the interpolation matrix for
% Semi-Lagrange for LFE.
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

nCoordinates =size(Mesh.Coordinates,1);

% Preallocate memory
pbE = zeros(3*nCoordinates,1);
I=zeros(3*nCoordinates,1);
J=zeros(3*nCoordinates,1);
level=[1 2 3];

% loop over all edges
for i = 1:nCoordinates
    
    p_vertex = tracedvertices(i,[1 2]);
    
    % Element for e1 and element mapping
    El_ID = tracedvertices(i,3);

    vid = Mesh.Elements(El_ID,:);
    a1 = Mesh.Coordinates(vid(1),:);
    a2 = Mesh.Coordinates(vid(2),:);
    a3 = Mesh.Coordinates(vid(3),:);

    % element mapping
    bK = a1;
    BK = [a2-bK; ...
        a3-bK];
    inv_BK = inv(BK);

    x_hat=(p_vertex-bK)*inv_BK;
   
    pbE(level(1)) = 1-sum(x_hat);
    pbE(level(2)) = x_hat(1);
    pbE(level(3)) = x_hat(2);
    I(level)=[i i i];
    J(level)=vid;
    level=level+3;

end% for
if(nargout > 1)
    varargout{1} = I;
    varargout{2} = J;
    varargout{3} = pbE;
else
    varargout{1} = sparse(I,J,pbE,nCoordinates,nCoordinates);
end
