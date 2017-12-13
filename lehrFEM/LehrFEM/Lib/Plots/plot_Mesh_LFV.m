function plot_Mesh_LFV(mesh)
	axes('DataAspectRatio',[1 1 1]);
	xl = min(mesh.Coordinates(:,1));
	xu = max(mesh.Coordinates(:,1));
	xp = .1*(xu-xl);
	yl = min(mesh.Coordinates(:,2));
	yu = max(mesh.Coordinates(:,2));
	yp = .1*(yu-yl);
	axis([xl-xp xu+xp yl-yp yu+yp]);
	hold on

	% Plot vertices
	for i=1:size(mesh.Coordinates,1)
		plot(mesh.Coordinates(i,1),mesh.Coordinates(i,2),'k*');
	end

	% Plot edges
	for i=1:size(mesh.Edges,1)
		plot([mesh.Coordinates(mesh.Edges(i,1),1) mesh.Coordinates(mesh.Edges(i,2),1)],...
		     [mesh.Coordinates(mesh.Edges(i,1),2) mesh.Coordinates(mesh.Edges(i,2),2)],'k--');
	end

	% Plot midpoints and centerpoints
	for i=1:size(mesh.Elements,1)
		mp1 = mesh.MidPoints(mesh.Vert2Edge(mesh.Elements(i,1),mesh.Elements(i,2)),:);
		mp2 = mesh.MidPoints(mesh.Vert2Edge(mesh.Elements(i,2),mesh.Elements(i,3)),:);
		mp3 = mesh.MidPoints(mesh.Vert2Edge(mesh.Elements(i,3),mesh.Elements(i,1)),:);
		cp = mesh.CenterPoints(i,:);

%		plot(mp1(1),mp1(2),'r');
%		plot(mp2(1),mp2(2),'r');
%		plot(mp3(1),mp3(2),'r');
%		plot(cp(1),cp(2),'r');
		plot([mp1(1) cp(1)],[mp1(2) cp(2)],'r-');
		plot([mp2(1) cp(1)],[mp2(2) cp(2)],'r-');
		plot([mp3(1) cp(1)],[mp3(2) cp(2)],'r-');
	end

end
