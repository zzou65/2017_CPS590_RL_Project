% MATLAB-LehrFEM script generating an unstructured triangular mesh for a
% circle domain using LehrFEM's built-in meshing tool.
% First define \Red{bounding box} for the 2D computational domain
% \Blue{$\Omega$} by specifying the lower left and upper right corners in
% the rows of a $2\times 2$-matrix.
bbox = [-1,-1;1,1];
% Largest admissible edge length of triangles; offers control for
% resolution of the mesh.
h0 = 0.3;
% Handle to \Red{signed distance function} \Blue{$\varphi(\Bx)$}:
% (distance from \Blue{$\partial\Omega$}, \Blue{$\varphi(\Bx)<0$} 
% $\Leftrightarrow$ \Blue{$\Bx\in\Omega$})  
dhd = @(x) sqrt(x(:,1).^2+x(:,2).^2)-1;
% Local element size function, here all set to $1$ $\Rightarrow$ request triangles
% of fairly uniform size. It allows to control the length of the edges
% locally.
hhandle = @(x) ones(size(x,1),1);
% Call the mesh generation function of LehrFEM.
Mesh = init_Mesh(bbox,h0,dhd,hhandle,[],1);
% Save the created mesh to file for later use
save_Mesh(Mesh,'coordinates.dat','elements.dat');
% Plot the mesh 
plot_Mesh(Mesh,'apts');
% Write the plot to an encapsulated PostScript file
print -depsc2 'circmeshex.eps';