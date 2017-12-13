function varargout = plotLine_LFE(U,Mesh,x_start,x_end)
% PLOTLINE_LFE Plot finite element solution.
%
%   PLOT_LFE(U,MESH,x_start,x_end) generates a plot of the finite element solution U on
%   the on the line x_start,x_end
%
%   The struct MESH must at least contain the following fields:
%    COORDINATES M-by-2 matrix specifying the vertices of the mesh.
%    ELEMENTS    N-by-3 matrix specifying the elements of the mesh.
%
%
%   Example:
%
%   plotLine_LFE(U,MESH,[0 0], [1 1]);

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
     nElements = size(Mesh.Elements,1);
     for i = 1:nElements
 
         % Extract vertices of current element
    
        idx = Mesh.Elements(i,:);
        a1 = Mesh.Coordinates(idx(1),:);
        a2 = Mesh.Coordinates(idx(2),:);
        a3 = Mesh.Coordinates(idx(3),:);
      
        % Compute element contributions
        s1_De=det([v;a2-a3]);
        if (s1_De~=0) 
            s1=det([v;a2-x])/s1_De;
            if (s1>=0 && s1<=1)
                ex=[ex,(a2+s1*(a3-a2)-x)*v'/norm(v)];
                eU=[eU,shap_LFE([1-s1,s1])*U(idx)];
            end
        end
        s2_De=det([v;a3-a1]);
        if (s2_De~=0) 
            s2=det([v;a3-x])/s2_De;
            if (s2>=0 && s2<=1)
                ex=[ex,(a3+s2*(a1-a3)-x)*v'/norm(v)];
                eU=[eU,shap_LFE([0,1-s2])*U(idx)];
            end
        end
        s3_De=det([v;a1-a2]);
        if (s3_De~=0) 
            s3=det([v;a1-x])/s3_De;
            if (s3>=0 && s3<=1)
                ex=[ex,(a1+s3*(a2-a1)-x)*v'/norm(v)];
                eU=[eU,shap_LFE([s3,0])*U(idx)];
            end
        end
     end       
     [ex,ID]=sort(ex);
     eU=eU(ID);
    
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