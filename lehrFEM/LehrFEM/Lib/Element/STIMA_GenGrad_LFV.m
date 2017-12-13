function aLoc = STIMA_GenGrad_LFV(vertices,midPoints,center,method,bdFlags,cHandle,kHandle,rHandle,conDom,varargin)
% STIMA_GENGRAD_LFV Element stiffness matrix for the generalized gradient.
%
%   ALOC = STIMA_GENLAPL_LFV(VERTICES,MIDPOINTS,CENTER,KHANDLE) computes the 
%   element stiffness matrix for the generalized gradient using linear
%   finite volumes.
%
%   VERTICES is 3-by-2 matrix specifying the vertices of the current element
%   in a row wise orientation.
%
%   MIDPOINTS is 3-by-2 matrix specifying the border points on the element
%   edges.
%
%   CENTER is a length 2 row vector specifying the center point of the element.
%
%   METHOD is 'barycentric' or 'orthogonal', the method used to make the mesh.
%
%   CHANDLE is a function handle for the c-function.
%
%   Example:
%
%   Aloc = STIMA_GenGrad_LFV([0 0; 1 0; 0 1],@(x)x);
%
%   Copyright 2007-2007 Eivind Fonn
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

	% Compute required coefficients

	m = zeros(3,3);
	m(1,2) = norm(midPoints(1,:)-center);
	m(1,3) = norm(midPoints(3,:)-center);
	m(2,3) = norm(midPoints(2,:)-center);
	m = m + m';

    %qr = gauleg(0, 1, 4, 1e-6);
    qr = gauleg(0, 1, 4);
    c = zeros(3,2);
    for i=1:length(qr.w)
        c(1,:) = c(1,:) + qr.w(i)*cHandle(midPoints(1,:)+qr.x(i)*(center-midPoints(1,:)));
        c(2,:) = c(2,:) + qr.w(i)*cHandle(midPoints(2,:)+qr.x(i)*(center-midPoints(2,:)));
        c(3,:) = c(3,:) + qr.w(i)*cHandle(midPoints(3,:)+qr.x(i)*(center-midPoints(3,:)));
    end

	gamma = zeros(3,3);
	gamma(1,2) = dot(confw(novec([midPoints(1,:); center]),vertices(2,:)-vertices(1,:)),c(1,:));
	gamma(2,3) = dot(confw(novec([midPoints(2,:); center]),vertices(3,:)-vertices(2,:)),c(2,:));
	gamma(1,3) = dot(confw(novec([midPoints(3,:); center]),vertices(3,:)-vertices(1,:)),c(3,:));
	gamma = gamma - gamma';

    if conDom
        mu = ones(3,3);
        mu(1,2) = kHandle((midPoints(1,:)+center)/2);
        mu(2,3) = kHandle((midPoints(2,:)+center)/2);
        mu(1,3) = kHandle((midPoints(3,:)+center)/2);
        mu = mu + mu';

        d = zeros(3,3);
        d(1,2) = norm(vertices(1,:)-vertices(2,:));
        d(2,3) = norm(vertices(2,:)-vertices(3,:));
        d(1,3) = norm(vertices(3,:)-vertices(1,:));
        d = d + d';

        z = gamma.*d./mu;
        r = zeros(3,3);
        for i=1:3, for j=1:3
            r(i,j) = rHandle(z(i,j));
        end; end;
    else
        r = ones(3,3)/2;
    end

	a = m.*gamma.*r;
	b = m.*gamma;

	% Compute stiffness matrix

	aLoc = b'-a';

	for i=1:3
		aLoc(i,i) = sum(a(i,:));
	end

	% Consider boundary edges

%	edges = [];
%	if bdFlags(1) < 0
%		edges = [edges; 1 2];
%	end
   
%	if bdFlags(2) < 0
%		edges = [edges; 2 3];
%	end
%	if bdFlags(3) < 0
%		edges = [edges; 1 3];
%	end
		
	% Iterate over boundary edges

%	for l=1:size(edges,1)
%		vids = edges(l,:);
%		i = vids(1);
%		j = vids(2);
%		vi = vertices(vids(1),:);
%		vj = vertices(vids(2),:);
%		if vids == [1 2]
%			ov = vertices(3,:);
%			mp = midPoints(1,:);
%		elseif vids == [2 3]
%			ov = vertices(1,:);
%			mp = midPoints(2,:);
%		else
%			ov = vertices(2,:);
%			mp = midPoints(3,:);
%		end
%
%		nu = confw(novec([vi; mp]),vi-ov); 
%		d = norm(vi-vj);

		% From i to j

%		eta = dot(nu,cHandle((vi+mp)/2));
%		s = norm(vi-mp);
%		aLoc(i,j) = aLoc(i,j) - (eta*s^2)/(2*d);
%		aLoc(i,i) = aLoc(i,i) - (eta*s^2)/(2*d);

		% From j to i

%		eta = dot(nu,cHandle((vj+mp)/2));
%		s = norm(vj-mp);
%		aLoc(j,i) = aLoc(j,i) - (eta*s^2)/(2*d);
%		aLoc(j,j) = aLoc(j,j) - (eta*s^2)/(2*d);
%	end

    aLoc = aLoc';

	% Helping functions

	function v = confw(v,t)
		if dot(v,t)<0
			v = -v;
		end
	end

	function out = novec(v)
		out = v(2,:)-v(1,:);
        if norm(out) ~= 0
            %norm(out)
            out = [-out(2) out(1)]/norm(out);
        end
	end
	
end
