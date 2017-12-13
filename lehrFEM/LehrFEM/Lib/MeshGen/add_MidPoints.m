function mesh = add_MidPoints(mesh,method)

    tol = 1e-6;

	% Input arguments

	if (nargin < 2) || (~strcmp(method,'barycentric') && ~strcmp(method,'orthogonal')) 
		method = 'barycentric';
	end

	% Get edge information and initialize

	if ~isfield(mesh,'Edges') || ~isfield(mesh,'Vert2Edge')
		mesh = add_Edges(mesh);
	end

	% Get center-points
	
	mesh.CenterPoints = zeros(size(mesh.Elements,1),2);
	for i=1:size(mesh.Elements,1)
		if strcmp(method, 'orthogonal')
			mesh.Type = 'orthogonal';
			a = (norm(mesh.Coordinates(mesh.Elements(i,2),:)-mesh.Coordinates(mesh.Elements(i,3),:)))^2;
			b = (norm(mesh.Coordinates(mesh.Elements(i,3),:)-mesh.Coordinates(mesh.Elements(i,1),:)))^2;
			c = (norm(mesh.Coordinates(mesh.Elements(i,1),:)-mesh.Coordinates(mesh.Elements(i,2),:)))^2;
			w(1) = a*(-a+b+c);
			w(2) = b*(a-b+c);
			w(3) = c*(a+b-c);
			if sum(w<-tol)>0
                w
				error('Obtuse triangle encountered.', 1);
			end
			w = w/sum(w); 
			mesh.CenterPoints(i,:) = w(1)*mesh.Coordinates(mesh.Elements(i,1),:)...
						+w(2)*mesh.Coordinates(mesh.Elements(i,2),:)...
						+w(3)*mesh.Coordinates(mesh.Elements(i,3),:);
		elseif strcmp(method, 'barycentric')
			mesh.Type = 'barycentric';
			mesh.CenterPoints(i,:) = sum(mesh.Coordinates(mesh.Elements(i,:),:),1)/size(mesh.Elements,2);
		end
	end

	% Get midpoints

	mesh.MidPoints = zeros(size(mesh.Edges,1),2);
	for i=1:size(mesh.Edges,1)
		mesh.MidPoints(i,:) = sum(mesh.Coordinates(mesh.Edges(i,:),:),1)/2;
	end

end
