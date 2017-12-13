function U = contour_DGCR(U,Mesh,varargin)
% contour_DGCR Contourplot for finite element solution.
%
%   CONTOUR_DGCR(U,MESH) generates a contourplot of the finite element
%   solution U on the mesh MESH
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 matrix specifying the elements of the mesh.
%
%   H = contour_DGCR(U,MESH) also returns the handle to the figure.
%
%   CONTOUR_LFE(U,MESH,LEVELS) draw LENGHT(LEVELS) contour lines at the
%   values specified in vector LEVELS.
%
%   CONTOUR_LFE(U,MESH,LEVELS,'c') lets you specify what color you want on
%   the level curves. See MATLAB help for colorspec for predefined colors.
%   Default color is blue.
%   
%   CONTOUR_LFE(U,MESH,LEVELS,'colorbar') appends a colorbar to the current
%   axes in the right location
%
%   Example:
%
%   contour_LFE(U,MESH);
%
%   See also contour

%   Copyright 2006-2006 Kari Borset
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

% Initialize constants


%Initialize constants

nCoordinates = size(Mesh.Coordinates,1);
nElements = size(Mesh.Elements,1);

xmin = min(Mesh.Coordinates(:,1));
xmax = max(Mesh.Coordinates(:,1));
xstart = xmin - 0.1*abs(xmin);
xend = xmax + 0.1*abs(xmax);
xlen = xend - xstart;

ymin = min(Mesh.Coordinates(:,2));
ymax = max(Mesh.Coordinates(:,2));
ystart = ymin - 0.1*abs(ymin);
yend = ymax + 0.1*abs(ymax);
ylen = yend - ystart;

hm = get_MeshWidth(Mesh);
len = min((xend-xstart),(yend-ystart));
h = get_MeshWidth(Mesh)/25;
h=max(h,len*sqrt(hm)/30);
hx = xlen/max(ceil(xlen/h));
hy = ylen/max(ceil(ylen/h));

%Generate grid

x=[xstart:hx:xend];
y=[ystart:hy:yend];
[X,Y]=meshgrid(x,y);
Z = zeros(size(X));

%Compute finiteelementssolution in each gridpoint

for i=1:nElements,
    
    %Extract vertex data
    
    Vertices = Mesh.Elements(i,:);
    P1 = Mesh.Coordinates(Vertices(1),:);
    P2 = Mesh.Coordinates(Vertices(2),:);
    P3 = Mesh.Coordinates(Vertices(3),:);
    
    %Create bounding box around element
    
    xmin = min([P1(1) P2(1) P3(1)]);
    xmax = max([P1(1) P2(1) P3(1)]);
    ymin = min([P1(2) P2(2) P3(2)]);
    ymax = max([P1(2) P2(2) P3(2)]);
    
    xi_s = floor((xmin-xstart)/hx)+1;
    xi_e = ceil((xmax-xstart)/hx)+1;
    yi_s = floor((ymin-ystart)/hy)+1;
    yi_e = ceil((ymax-ystart)/hy)+1;
    for k = xi_s:xi_e,
        for l = yi_s:yi_e,
        
            %Compute Element mapping
            
            bK = P1;
            BK = [P2-bK; P3-bK];
            inv_BK = 1/(BK(1,1)*BK(2,2)-BK(1,2)*BK(2,1))*[BK(2,2) -BK(1,2) ; -BK(2,1) BK(1,1)];
            
            %Transform grid points back to reference element
            
            Px(1) = X(l,k);
            Px(2) = Y(l,k);
            xhat = (Px-bK)*inv_BK;
            
            %Check if gridpoint is inside element
            
            if(xhat(1) >= 0 && xhat(2) >= 0 && xhat(1)+xhat(2) <= 1),
                %Compute solution in gridpoint
                N = shap_DGCR(xhat);
                Z(l,k) = U(i*3-2)*N(1) + ...
                         U(i*3-1)*N(2) + ...
                         U(i*3)*N(3);
            end                     
        end
    end
end

if nargin == 3,
   levels = varargin{1};
   color = 'b';
   
   fig = figure('Name','Crouzeix-Raviart finite elements');
   contour(X,Y,Z,levels,color);
   grid('on');
   set(gca,'DataAspectRatio',[1 1 1]);
   
elseif nargin == 4,
    levels = varargin{1};
    if varargin{2} == 'colorbar',
        fig = figure('Name','Crouzeix-Raviart finite elements');
        contour(X,Y,Z,levels);
        colorbar;
        grid('on');
        set(gca,'DataAspectRatio',[1 1 1]);
    else
        color = varargin{2};

        fig = figure('Name','Crouzeix-Raviart finite elements');
        contour(X,Y,Z,levels,color);
        grid('on');
        set(gca,'DataAspectRatio',[1 1 1]);
    end
    
else
    color = 'b';
    fig = figure('Name','Linear finite elements');
    contour(X,Y,Z,color);
    grid('on');
    set(gca,'DataAspectRatio',[1 1 1]);
   
end

if(nargout > 0)
  varargout{1} = fig;  
end
    
return


