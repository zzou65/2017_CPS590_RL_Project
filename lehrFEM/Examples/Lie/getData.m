function d=getData(n);

switch n
    %
    case 1
        % velocity v
        data.V_Handle=@(x,varargin)[ones(size(x,1),1) 0.5.*ones(size(x,1),1)];
        % solution u
        data.U_EX_Handle=@(x,varargin)[(1-x(:,2)).*(1+x(:,2)).*sin(pi*x(:,1)) ...
                                (1-x(:,1)).*(1+x(:,1)).*sin(pi*x(:,2))];
        %  v *u                  
        data.V_U_EX_Handle=@(x,varargin)(1-x(:,2)).*(1+x(:,2)).*sin(pi*x(:,1))+...
                                    0.5*(1-x(:,1)).*(1+x(:,1)).*sin(pi*x(:,2)) ;
        %  - v x curl u                          
        data.SOL1_Handle=@(x,varargin)[  -x(:,2).*sin(pi.*x(:,1))+   x(:,1).*sin(pi.*x(:,2)) ...
                                      +2.*x(:,2).*sin(pi.*x(:,1))-2.*x(:,1).*sin(pi.*x(:,2))];    % - v x curl u
        % grad(v*u)
        data.SOL2_Handle=@(x,varargin)[-pi.*(-1+x(:,2).^2).*cos(pi.*x(:,1))-x(:,1).*sin(pi.*x(:,2)) ...
                                       -pi/2.*(-1+x(:,1).^2).*cos(pi.*x(:,2))-2*x(:,2).*sin(pi.*x(:,1))]; % grad (v u)
        % curl curl u
        data.SOL3_Handle=@(x,varargin)[2.*(-pi.*x(:,1).*cos(pi.*x(:,2))+sin(pi.*x(:,1))) ...
                                       2.*(-pi.*x(:,2).*cos(pi.*x(:,1))+sin(pi.*x(:,2)))]; % curl curl u
        % curl u
        data.CURL_U_EX_Handle=@(x,varargin)-2.*x(:,2).*sin(pi.*x(:,1))+2.*x(:,1).*sin(pi.*x(:,2));    % curl u
        
        data.boundtype=[-1 -1 -2 -2];
        data.Coordinates=[-1 -1; 1 -1; 1 1; -1 1];
        data.Elements=[1 2 4;2 3 4];
        
     case 2
        % velocity v
        data.V_Handle=@(x,varargin)ones(size(x,1),2);
        % solution u
        data.U_EX_Handle=@(x,varargin)[sin(pi.*x(:,2)).*sin(pi.*x(:,1)) ...
                                       sin(pi.*x(:,1)).*sin(pi.*x(:,2))];
        %  v *u                  
        data.V_U_EX_Handle=@(x,varargin)(sin(pi.*x(:,2)).*sin(pi.*x(:,1))+sin(pi.*x(:,1)).*sin(pi.*x(:,2)));
        %  - v x curl u                          
        data.SOL1_Handle=@(x,varargin)[+pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1))-pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2)) ...
                                       -pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1))+pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2))]; 
        % grad(v*u)
        data.SOL2_Handle=@(x,varargin)[2.*pi.*sin(pi.*x(:,2)).*cos(pi.*x(:,1))...
                                       2.*pi.*sin(pi.*x(:,1)).*cos(pi.*x(:,2))]; 
        % curl curl u
        data.SOL3_Handle=@(x,varargin)[pi.^2.*cos(pi.*x(:,2)).*cos(pi.*x(:,1))+pi.^2.*sin(pi.*x(:,1)).*sin(pi.*x(:,2)) ...
                                      +pi.^2.*cos(pi.*x(:,2)).*cos(pi.*x(:,1))+pi.^2.*sin(pi.*x(:,1)).*sin(pi.*x(:,2))];
        % curl u
        data.CURL_U_EX_Handle=@(x,varargin)(pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1))-pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2)));
        
        data.boundtype=[-1 -1 -2 -2];
        data.Coordinates=[-1 -1; 1 -1; 1 1; -1 1];
        data.Elements=[1 2 4;2 3 4];
     
    case 3
        % velocity v
        data.V_Handle=@(x,varargin)ones(size(x,1),2);
        % solution u
        data.U_EX_Handle=@(x,varargin)[1-x(:,1).^2 1-x(:,2).^2];
        %  v *u                  
        data.V_U_EX_Handle=@(x,varargin)2-x(:,1).^2-x(:,2).^2;
        %  - v x curl u                          
        data.SOL1_Handle=@(x,varargin)zeros(size(x,1),2);
        % grad(v*u)
        data.SOL2_Handle=@(x,varargin)[-2.*x(:,1) -2.*x(:,2)];
        % curl curl u
        data.SOL3_Handle=@(x,varargin)zeros(size(x,1), 2);
        % curl u
        data.CURL_U_EX_Handle=@(x,varargin)zeros(size(x,1),1);            
        data.boundtype=[-1 -1 -2 -2];
        data.Coordinates=[-1 -1; 1 -1; 1 1; -1 1];
        data.Elements=[1 2 4;2 3 4];  
    
    case 4
        % velocity v
        data.V_Handle=@(x,varargin)[ones(size(x,1),1) ones(size(x,1),1)];
        % solution u
        data.U_EX_Handle=@(x,varargin)[1-x(:,1) 1-x(:,2)];
        %  v *u                  
        data.V_U_EX_Handle=@(x,varargin)2-x(:,1)-x(:,2);
        %  - v x curl u                          
        data.SOL1_Handle=@(x,varargin)zeros(size(x,1),2);
        % grad(v*u)
        data.SOL2_Handle=@(x,varargin)-ones(size(x(:,1),1),2);
        % curl curl u
        data.SOL3_Handle=@(x,varargin)zeros(size(x,1),2)  ;
        % curl u
        data.CURL_U_EX_Handle=@(x,varargin)zeros(size(x,1),1) ;
        data.boundtype=[-1 -1 -2 -2];
        data.Coordinates=[-1 -1; 1 -1; 1 1; -1 1];
        data.Elements=[1 2 4;2 3 4];
    
    case 5
        % velocity v
        data.V_Handle=@(x,varargin)[(2+sin(pi*x(:,1))).*ones(size(x,1),1) (2+sin(pi*x(:,1))).*ones(size(x,1),1)];
        % solution u
        data.U_EX_Handle=@(x,varargin)[(1-x(:,2)).*(1+x(:,2)).*sin(pi*x(:,1)) ...
                                       (1-x(:,1)).*(1+x(:,1)).*sin(pi*x(:,2))];
        %  v *u                  
        data.V_U_EX_Handle=@(x,varargin)[(1-x(:,2).^2).*sin(pi*x(:,1)).*(2+sin(pi*x(:,1)))+...
                                         (1-x(:,1).^2).*sin(pi*x(:,2)).*(2+sin(pi*x(:,1)))];
        %  - v x curl u                          
        data.SOL1_Handle=@(x,varargin)[(2.*(2+sin(pi*x(:,1))).*(-x(:,2).*sin(pi*x(:,1))+x(:,1).*sin(pi*x(:,2)))) ...
                                       (2.*(2+sin(pi*x(:,1))).*(+x(:,2).*sin(pi*x(:,1))-x(:,1).*sin(pi*x(:,2))))];
        % grad(v*u)
        data.SOL2_Handle=@(x,varargin)[-pi*(-1+x(:,2).^2).*(2.*cos(pi.*x(:,1))+sin(2.*pi.*x(:,1)))-(pi.*(-1+x(:,1).^2).*cos(pi.*x(:,1))+2.*x(:,1).*(2+sin(pi.*x(:,1)))).*sin(pi.*x(:,2)) ...
                                       (2+sin(pi.*x(:,1))).*(-pi.*(-1+x(:,1).^2).*cos(pi.*x(:,2))-2.*x(:,2).*sin(pi.*x(:,1)))]
        % curl curl u
        data.SOL3_Handle=@(x,varargin)[2.*(-pi.*x(:,1).*cos(pi.*x(:,2))+sin(pi.*x(:,1))) ...
                                       2.*(-pi.*x(:,2).*cos(pi.*x(:,1))+sin(pi.*x(:,2)))];
        % curl u
        data.CURL_U_EX_Handle=@(x,varargin)-2.*x(:,2).*sin(pi.*x(:,1))+2.*x(:,1).*sin(pi.*x(:,2));
        
        data.boundtype=[-1 -1 -2 -2];
        data.Coordinates=[-1 -1; 1 -1; 1 1; -1 1];
        data.Elements=[1 2 4;2 3 4];
        
    case 6
        % velocity v
        data.V_Handle=@(x,varargin)[-ones(size(x,1),1)-x(:,2) ones(size(x,1),1)+x(:,1)];
        % solution u
        data.U_EX_Handle=@(x,varargin)[(1-x(:,2)).*(1+x(:,2)).*sin(pi*x(:,1)) ...
                                       (1-x(:,1)).*(1+x(:,1)).*sin(pi*x(:,2))];
        %  v *u                  
        data.V_U_EX_Handle=@(x,varargin)(-1-x(:,2)).*(1-x(:,2).^2).*sin(pi*x(:,1))+...
                                        (1+x(:,1)).*(1-x(:,1).^2).*sin(pi*x(:,2)) ;
        %  - v x curl u                          
        data.SOL1_Handle=@(x,varargin)[(2.*(1+x(:,1)).*(-x(:,2).*sin(pi.*x(:,1))+x(:,1).*sin(pi.*x(:,2)))) ...
                                       (2.*(1+x(:,2)).*(-x(:,2).*sin(pi.*x(:,1))+x(:,1).*sin(pi.*x(:,2))))];
        % grad(v*u)
        data.SOL2_Handle=@(x,varargin)[-pi.*(1-x(:,2).^2).*(1+x(:,2)).*cos(pi.*x(:,1))+(1-2*x(:,1)-3.*x(:,1).^2).*sin(pi.*x(:,2)) ...
                                       +pi.*(1-x(:,1).^2).*(1+x(:,1)).*cos(pi.*x(:,2))-(1-2*x(:,2)-3.*x(:,2).^2).*sin(pi.*x(:,1))];
        % curl curl u
        data.SOL3_Handle=@(x,varargin)[2.*(-pi.*x(:,1).*cos(pi.*x(:,2))+sin(pi.*x(:,1))) ...
                                       2.*(-pi.*x(:,2).*cos(pi.*x(:,1))+sin(pi.*x(:,2)))];
        % curl u
        data.CURL_U_EX_Handle=@(x,varargin)2.*x(:,2).*sin(pi.*x(:,1))-2.*x(:,1).*sin(pi.*x(:,2));
        
        data.boundtype=[-1 -2 -1 -2];
        data.Coordinates=[-1 -1; 1 -1; 1 1; -1 1];
        data.Elements=[1 2 4;2 3 4];
end
        
        
d=data;
return