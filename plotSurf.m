function Fhandel=plotSurf(points, value, crange, cbar_on)
x=points(:,1);
y=points(:,2);
z=value;
%plot3(x,y,z,'.-')
tri = delaunay(x,y);
Fhandel = trisurf(tri, x, y, z);
%axis vis3d
axis on
axis equal
%lighting phong
shading interp
if(~isempty(crange))
    caxis(crange)
end
if(cbar_on==1)
    h=colorbar;
    set(h, 'fontsize', 20)
end
view([0,0,90])
