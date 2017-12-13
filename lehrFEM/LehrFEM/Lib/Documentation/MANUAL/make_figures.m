% Generates .eps files for the manual
%
%   Possible values for FLAG are:
%    0 Generate .eps files of mesh plotting routines.
%    1 Generate .eps files for uniform red refinements.
%    2 Genearte .eps files for largest edge bisection.
%    3 Generate .eps files for example of largest edge bisection

switch(flag)
  case 1  
    
    % .eps files of mesh plotting routines

      % Initialize constants
      
      BBox = [0 0; 1 1];
      h0 = 0.05;
      DHandle = inline('dist_diff(dist_diff(dist_diff(dist_circ(x,[0 0],1),x(:,2)),x(:,1)),dist_circ(x,[0 0],.5))','x');
      HHandle = @h_uniform;
      FixedPos = [0 .5; .5 0; 1 0; 0 1];
      disp = 1;  
      EPS_1 = 'mesh.eps';
      EPS_2 = 'qual.eps';
      EPS_3 = 'usr.eps';
  
      % Generate mesh
  
      Mesh = init_Mesh(BBox,h0,DHandle,HHandle,FixedPos,disp);
  
      % Generate plots of the mesh
  
      plot_Mesh(Mesh);
      print('-deps',EPS_1);
      system(['gv ',EPS_1,' &']);
  
      plot_Qual(Mesh);
      print('-depsc',EPS_2)
      system(['gv ',EPS_2,' &']);
  
      plot_USR(Mesh);
      print('-deps',EPS_3);
      system(['gv ',EPS_3,' &']);
  
      % Clear memory and close all leftover figures
      
      close all;
      clear all;

  case 2

    % Generates the .eps files for regular red refinements

      % Initialize constants
  
      SHIFT_1 = 0.1;
      SHIFT_2 = 0.2;
      FONT_SIZE = 16;
      FONT_WEIGHT = 'bold';
      FONT_COLOR = 'k';
      EPS_1 = 'REG_org.eps';
      EPS_2 = 'REG_ref.eps';
  
      % Vertices 

      Coordinates = [1/2 1/2; ...
                       2 1/4; ...
                       1   2]; 
      Coordinates(4,:) = 1/2*(Coordinates(1,:)+Coordinates(2,:));
      Coordinates(5,:) = 1/2*(Coordinates(2,:)+Coordinates(3,:));
      Coordinates(6,:) = 1/2*(Coordinates(3,:)+Coordinates(1,:));
  
      % Definition of elements
  
      Elem_Old = [1 2 3];
      Elem_New = [1 4 6; ...
                  4 2 5; ...
                  6 5 3; ...
                  4 5 6];
  
      % Compute midpoints of elements
  
      Mid_Old = 1/3*(Coordinates(Elem_Old(1),:)+Coordinates(Elem_Old(2),:)+Coordinates(Elem_Old(3),:));
      Mid_New = [1/3*(Coordinates(Elem_New(1,1),:)+Coordinates(Elem_New(1,2),:)+Coordinates(Elem_New(1,3),:)); ...
                 1/3*(Coordinates(Elem_New(2,1),:)+Coordinates(Elem_New(2,2),:)+Coordinates(Elem_New(2,3),:)); ...
                 1/3*(Coordinates(Elem_New(3,1),:)+Coordinates(Elem_New(3,2),:)+Coordinates(Elem_New(3,3),:))
                 1/3*(Coordinates(Elem_New(4,1),:)+Coordinates(Elem_New(4,2),:)+Coordinates(Elem_New(4,3),:))];
  
      % Generate figures
  
      h1 = figure;
      patch('Faces',Elem_Old, ...
            'Vertices',Coordinates, ...
            'FaceColor','none', ...
            'EdgeColor','k');
      for i = 1:3
        Shifted_Coord = Coordinates(Elem_Old(i),:) + SHIFT_1*(Mid_Old-Coordinates(Elem_Old(i),:));
        text(Shifted_Coord(1),Shifted_Coord(2), int2str(i), ...
             'HorizontalAlignment', 'Center', ...
             'VerticalAlignment', 'Middle',...
             'FontSize', FONT_SIZE, ...
             'FontWeight',FONT_WEIGHT, ...
             'Color', FONT_COLOR);
      end
      axis('equal');
      axis('off');
      print('-deps',EPS_1);
      close(h1); 
      system(['gv ',EPS_1,' &']);  
  
      h2 = figure;
      patch('faces',Elem_New, ...
            'vertices',Coordinates, ...
            'facecolor','none', ...
            'edgecolor','k');
      for j = 1:4
        for i = 1:3
          Shifted_Coord = Coordinates(Elem_New(j,i),:) + SHIFT_2*(Mid_New(j,:)-Coordinates(Elem_New(j,i),:));
          text(Shifted_Coord(1),Shifted_Coord(2), int2str(i), ...
             'HorizontalAlignment', 'Center', ...
             'VerticalAlignment', 'Middle',...
             'FontSize', FONT_SIZE, ...
             'FontWeight',FONT_WEIGHT, ...
             'Color', FONT_COLOR);
        end
      end
      axis('equal');
      axis('off');
      print('-deps',EPS_2);
      close(h2);
      system(['gv ',EPS_2,' &']);
  
      % Clear memory
  
      clear all;

  case 3

    % Generates the .eps files for largest edge bisection

      % Initialize constants
 
      SHIFT_1 = 0.1;
      SHIFT_2 = 0.2;
      FONT_SIZE = 16;
      FONT_WEIGHT = 'bold';
      FONT_COLOR = 'k';
      EPS_1 = 'LEB_org.eps';
      EPS_2 = 'LEB_ref.eps';
  
      % Vertices 

      Coordinates = [1/2 1/2; ...
                       2 1/4; ...
                       1   2; ...
                       2 3/2]; 
      Coordinates(5,:) = 1/2*(Coordinates(2,:)+Coordinates(3,:));
  
      % Definition of elements
  
      Elem_Old = [2 3 1; ... 
                  3 2 4];
      Elem_New = [1 2 5; ...
                  2 4 5; ...
                  4 3 5; ...
                  3 1 5];
  
      % Compute midpoints of elements
  
      Mid_Old = [1/3*(Coordinates(Elem_Old(1,1),:)+Coordinates(Elem_Old(1,2),:)+Coordinates(Elem_Old(1,3),:)); ...
                 1/3*(Coordinates(Elem_Old(2,1),:)+Coordinates(Elem_Old(2,2),:)+Coordinates(Elem_Old(2,3),:))];
      Mid_New = [1/3*(Coordinates(Elem_New(1,1),:)+Coordinates(Elem_New(1,2),:)+Coordinates(Elem_New(1,3),:)); ...
                 1/3*(Coordinates(Elem_New(2,1),:)+Coordinates(Elem_New(2,2),:)+Coordinates(Elem_New(2,3),:)); ...
                 1/3*(Coordinates(Elem_New(3,1),:)+Coordinates(Elem_New(3,2),:)+Coordinates(Elem_New(3,3),:))
                 1/3*(Coordinates(Elem_New(4,1),:)+Coordinates(Elem_New(4,2),:)+Coordinates(Elem_New(4,3),:))];
         
      % Generate figures
  
      h1 = figure;
      patch('Faces',Elem_Old, ...
            'Vertices',Coordinates, ...
            'FaceColor','none', ...
            'EdgeColor','k');
      for j = 1:2  
        for i = 1:3
          Shifted_Coord = Coordinates(Elem_Old(j,i),:) + SHIFT_1*(Mid_Old(j,:)-Coordinates(Elem_Old(j,i),:));
          text(Shifted_Coord(1),Shifted_Coord(2), int2str(i), ...
               'HorizontalAlignment', 'Center', ...
               'VerticalAlignment', 'Middle',...
               'FontSize', FONT_SIZE, ...
               'FontWeight',FONT_WEIGHT, ...
               'Color', FONT_COLOR);
        end
      end
      axis('equal');
      axis('off');
      print('-deps',EPS_1);
      close(h1);
      system(['gv ',EPS_1,' &']);  
  
      h2 = figure;
      patch('faces',Elem_New, ...
            'vertices',Coordinates, ...
            'facecolor','none', ...
            'edgecolor','k');
      for j = 1:4
        for i = 1:3
          Shifted_Coord = Coordinates(Elem_New(j,i),:) + SHIFT_2*(Mid_New(j,:)-Coordinates(Elem_New(j,i),:));
          text(Shifted_Coord(1),Shifted_Coord(2), int2str(i), ...
             'HorizontalAlignment', 'Center', ...
             'VerticalAlignment', 'Middle',...
             'FontSize', FONT_SIZE, ...
             'FontWeight',FONT_WEIGHT, ...
             'Color', FONT_COLOR);
        end
      end
      axis('equal');
      axis('off');
      print('-deps',EPS_2);
      close(h2);
      system(['gv ',EPS_2,' &']);
  
      % Clear memory
  
      clear all;
  
      
  case 4

    % Initialize constants    
      
    NREFS = 15;
    EPS_1 = 'LShap.eps';
    EPS_2 = 'LShap_LEB.eps';
  
    % Initialze mesh
    
    Mesh.Coordinates = [-1 -1; ...
                         0 -1; ...
                         0  0; ...
                         1  0; ...
                         1  1; ...
                        -1  1];
    Mesh.Elements = [1 2 3; ...
                     3 4 5;
                     3 5 6; ...
                     1 3 6];
    plot_Mesh(Mesh);
                     
    % Add edge data structure
  
    Mesh = add_Edges(Mesh);
    Loc = get_BdEdges(Mesh);
    Mesh.BdFlags = zeros(size(Mesh.Edges,1),1);
    Mesh.BdFlags(Loc) = -1;
 
    % Initialize largest edge bisection
  
    Mesh = init_LEB(Mesh);  
    print('-deps',EPS_1);
    
    system(['gv ',EPS_1,' &']);
  
    % Do largest edge bisection
  
    for i = 1:NREFS
      Mesh = refine_LEB(Mesh,1);
    end
    plot_Mesh(Mesh); 
    print('-deps',EPS_2);
    
    system(['gv ',EPS_2,' &']);
    
    % Clear memory and close left over figures
    
    close all;
    clear all;
end