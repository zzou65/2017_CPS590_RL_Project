function varargout = plotLine_QFE(U,Mesh,x_start,x_end)
% PLOTLINE_QFE Plot finite element solution.
%
%   PLOT_QFE(U,MESH,x_start,x_end) generates a plot of the finite element solution U on
%   the on the line x_start,x_end
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 matrix specifying the elements of the mesh.
%
%
%   Example:
%
%   plotLine_QFE(U,MESH,[0 0], [1 1]);

%   Copyright 2008-2008 Holger Heumann
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

  % Initialize constants
  
  % Generate figure
     ex=[];
     eU=[];
     x=x_start;
     v=x_end-x_start;
     nCoordinates = size(Mesh.Coordinates,1);
     nElements = size(Mesh.Elements,1);
     for i = 1:nElements
 
         % Extract vertices of current element
        l=[];
        vidx = Mesh.Elements(i,:);    
        idx = [vidx ...
           Mesh.Vert2Edge(Mesh.Elements(i,1),Mesh.Elements(i,2))+nCoordinates ...
           Mesh.Vert2Edge(Mesh.Elements(i,2),Mesh.Elements(i,3))+nCoordinates ...
           Mesh.Vert2Edge(Mesh.Elements(i,3),Mesh.Elements(i,1))+nCoordinates];
        
        a1 = Mesh.Coordinates(idx(1),:);
        a2 = Mesh.Coordinates(idx(2),:);
        a3 = Mesh.Coordinates(idx(3),:);
      
        % Compute element contributions
        s1_De=det([v;a2-a3]);
        if (s1_De~=0) 
            s1=det([v;a2-x])/s1_De;
            if (s1>=0 && s1<=1)
                ex=[ex;(a2+s1*(a3-a2)-x)*v'/norm(v)];
                eU=[eU;shap_QFE([1-s1 s1])*U(idx,:)];
                l=[l;1-s1 s1];
            end
        end
        s2_De=det([v;a3-a1]);
        if (s2_De~=0) 
            s2=det([v;a3-x])/s2_De;
            if (s2>=0 && s2<=1)
                ex=[ex;(a3+s2*(a1-a3)-x)*v'/norm(v)];
                eU=[eU;shap_QFE([0 1-s2])*U(idx,:)];
                l=[l;0 1-s2];
            end
        end
        s3_De=det([v;a1-a2]);
        if (s3_De~=0) 
            s3=det([v;a1-x])/s3_De;
            if (s3>=0 && s3<=1)
                ex=[ex;(a1+s3*(a2-a1)-x)*v'/norm(v)];
                eU=[eU;shap_QFE([s3 0])*U(idx,:)];
                l=[l ;s3 0];
            end
        end
        if (size(l,1)==2) 
            ex=[ex;ex(end-1)+1/2*(ex(end)-ex(end-1))];
            eU=[eU;shap_QFE(l(1,:)+1/2*(l(2,:)-l(1,:)))*U(idx,:)];
        end
            
     end       
     
     [ex,ID]=sort(ex);
     eU=eU(ID,:);
     m=1;  
     for i=1:(max(size(ID))-1);
       if (abs(ex(i)-ex(i+1))>eps)
           m=[m,i+1];
       end
     end
     %[ex,m,n]=unique(ex);
     ex=ex(m);
     eU=eU(m,:);
     
    
     %Compute axes limits
  
%      YMin =min(eU);
%      YMax =max(eU);
%      if(YMin < YMax)
%        YLim = [YMin YMax] + OFFSET*(YMax-YMin)*[-1 1];
%      else
%       YLim = [1-OFFSET 1+OFFSET]*YMin;   
%      end

     % Generate figure
  
    if(nargout > 0)
      varargout{1} = ex;
      varargout{2} = eU;
    else
     fig = figure('Name','profile plot');
     plot(ex,eU,'r-');
    end
return