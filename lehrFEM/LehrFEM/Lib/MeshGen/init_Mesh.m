function Mesh = init_Mesh(BBox,h0,DHandle,HHandle,FixedPos,disp,varargin)
% INIT_MESH Mesh generator.
%
%   MESH = INIT_MESH(BBOX,H0,DHANDLE,HHANDLE,DISP,FIXEDPOS) generates a triangular
%   mesh of the geometry specified by the signed distance function DHANDLE (function
%   handle/inline object).
%   
%   The bounding box BBOX and the initial mesh width H0 are used to produce
%   an initial triangulation of the domain. The first and second column of
%   BBOX serve as axes limits for the x and y axes. 
%
%   The function HHANDLE (function handle/inline object) specifies the desired
%   relative element size distribution. The element size distribution function
%   should not be zero anywhere inside the domain.
%   
%   FIXEDPOS specifies the number of fixed points on the boundary of the
%   domain. When using heuristics of signed distance functions it is
%   important to fix the location of the corner points of the domain in
%   order to obtain good results
%
%   The display flag DISP specifies wheter progress should be displayed
%   during the meshing procedure:
%    1 The mesh is redrawn after every delaunay triangulation
%    0 The mesh will not be displayed during the meshing procedure
%
%   MESH = GEN_MESH(BBOX,H0,DHANDLE,HHANDLE,FIXEDPOS,DISP,FPARAM) also handles
%   the variable argument list FPARAM to the functions DHANDLE and HHANDLE.
%
%   Example:
%
%   BBox = [-3 -3; 3 3], h0 = 0.3, FixedPos = []
%   DHandle = @dist_circ, C = [0 0], R = 3
%   HHandle = @h_uniform
%   disp = 1
%   Mesh = gen_Mesh(BBox,h0,DHandle,HHandle,FixedPos,disp,C,R);

%   Copyright 2004-2005 Peter-Olof Persson
%   Department of Mathematic
%   Massachusetts Institute of Technology
%   Massachusetts Avenue 77
%   Cambridge MA 02139

  % Initialize constants

  DPTOL = 0.001;        % Tolerance for stopping criterion
  TTOL = 0.1;           % Retriangulation tolerance
  FSCALE = 1.2;         % Force scaling factor
  DELTAT = 0.2;         % Time step size
  GEPS = 0.001*h0;      % Geometric resolution
  DEPS = sqrt(eps)*h0;  % Minimum Geometric resolution
  OFFSET = 0.05;        % Offset parameter
  
  % Compute axis limits and open figure
  
  if(disp)
    fig = figure('Name','Unstructured mesh generator');
    XLim = [BBox(1,1) BBox(2,1)] + (BBox(2,1)-BBox(1,1))*OFFSET*[-1 1];
    YLim = [BBox(1,2) BBox(2,2)] + (BBox(2,2)-BBox(1,2))*OFFSET*[-1 1];
  end
  
  % Create initial distribution in bounding box
  
  [x,y] = meshgrid(BBox(1,1):h0:BBox(2,1),BBox(1,2):h0*sqrt(3)/2:BBox(2,2));
  x(2:2:end,:) = x(2:2:end,:)+h0/2;
  p = [x(:),y(:)];

  % Remove points outside the region, apply the rejection method
  
  p = p(feval(DHandle,p,varargin{:}) < GEPS,:);
  r0 = 1./feval(HHandle,p,varargin{:}).^2;
  p = p(rand(size(p,1),1) < r0./max(r0),:);
  for i = 1:size(FixedPos,1)
    p = p(sqrt((p(:,1)-FixedPos(i,1)).^2+(p(:,2)-FixedPos(i,2)).^2) >= h0,:);
  end
  p = [FixedPos; p];
  N = size(p,1);
  
  pold=inf;
  while(1)
      
    % Retriangulate using a delaunay tessellation
      
    if(max(sqrt(sum((p-pold).^2,2))/h0) > TTOL)
      pold = p;
      t = delaunayn(p);                                  
      pmid = (p(t(:,1),:)+p(t(:,2),:)+p(t(:,3),:))/3;
      t = t(feval(DHandle,pmid,varargin{:}) < -GEPS,:);
      
      % Describe each bar by a unique pair of nodes
      
      bars = [t(:,[1,2]); t(:,[1,3]); t(:,[2,3])];
      bars = unique(sort(bars,2),'rows'); 
      
      % Display meshing progress
      
      if(disp)
        clf;
        patch('faces',t, ...
              'vertices',[p(:,1) p(:,2)], ...
              'facecolor','none', ...
              'edgecolor','k');
        set(gca,'XLim',XLim,'YLim',YLim,'DataAspectRatio',[1 1 1],'Box','on');
        drawnow;
      end
    end

    % Move mesh points based on bar lengths L and forces F
    
    barvec = p(bars(:,1),:)-p(bars(:,2),:);
    L = sqrt(sum(barvec.^2,2));
    hbars = feval(HHandle,(p(bars(:,1),:)+p(bars(:,2),:))/2,varargin{:});
    L0 = hbars*FSCALE*sqrt(sum(L.^2)/sum(hbars.^2));
    F = max(L0-L,0);
    Fvec = F./L*[1,1].*barvec;
    Ftot = full(sparse(bars(:,[1,1,2,2]),ones(size(F))*[1,2,1,2],[Fvec,-Fvec],N,2));
    Ftot(1:size(FixedPos,1),:)=0;
    p = p+DELTAT*Ftot;

    % Bring outside points back to the boundary
    
    d = feval(DHandle,p,varargin{:});
    ix = d > 0;
    dgradx = (feval(DHandle,[p(ix,1)+DEPS,p(ix,2)],varargin{:})-d(ix))/DEPS;
    dgrady = (feval(DHandle,[p(ix,1),p(ix,2)+DEPS],varargin{:})-d(ix))/DEPS;
    p(ix,:) = p(ix,:)-[d(ix).*dgradx,d(ix).*dgrady];

    % Termination criterion: All interior nodes move less than DPTOL
    
    dp = max(sqrt(sum(DELTAT*Ftot(d<-GEPS,:).^2,2))/h0);
    if(dp < DPTOL)
      break;
    end
  end  
  
  % Close all open figures
  
  if(disp)
    close(fig);
  end
  
  % Assign output arguments
  
  Mesh.Coordinates = p;
  Mesh.Elements = t;
  
  % Change orientation of all elements to counter-clockwise
  
  for i = 1:size(Mesh.Elements,1);
    v = Mesh.Coordinates(Mesh.Elements(i,2),:)-Mesh.Coordinates(Mesh.Elements(i,1),:);
    w = Mesh.Coordinates(Mesh.Elements(i,3),:)-Mesh.Coordinates(Mesh.Elements(i,1),:);
    if(v(1)*w(2)-v(2)*w(1) < 0)
      tmp = Mesh.Elements(i,1);
      Mesh.Elements(i,1) = Mesh.Elements(i,2);
      Mesh.Elements(i,2) = tmp;
    end
  end
  
return
