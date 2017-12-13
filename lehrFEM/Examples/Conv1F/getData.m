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
                                    0.5*(1-x(:,1)).*(1+x(:,1)).*sin(pi*x(:,2));
        %  - v x curl u                          
        data.SOL1_Handle=@(x,varargin)[  -x(:,2).*sin(pi.*x(:,1))+x(:,1).*sin(pi.*x(:,2)) ...
                                      +2.*x(:,2).*sin(pi.*x(:,1))-2.*x(:,1).*sin(pi.*x(:,2))];    % - v x curl u
        % grad(v*u)
        data.SOL2_Handle=@(x,varargin)[-pi.*(-1+x(:,2).^2).*cos(pi.*x(:,1))-x(:,1).*sin(pi.*x(:,2)) ...
                                       -pi/2.*(-1+x(:,1).^2).*cos(pi.*x(:,2))-2*x(:,2).*sin(pi.*x(:,1))]; % grad (v u)
        % curl curl u
        data.SOL3_Handle=@(x,varargin)[2.*(-pi.*x(:,1).*cos(pi.*x(:,2))+sin(pi.*x(:,1))) ...
                                       2.*(-pi.*x(:,2).*cos(pi.*x(:,1))+sin(pi.*x(:,2)))]; % curl curl u
        % - curl u
        data.CURL_U_EX_Handle=@(x,varargin) -2.*x(:,2).*sin(pi.*x(:,1))+2.*x(:,1).*sin(pi.*x(:,2));    % curl u
        
        % div u
        data.Div_U_EX_Handle=@(x,varargin)[pi*(1-x(:,2)).*(1+x(:,2)).*cos(pi*x(:,1))+...
                                pi*(1-x(:,1)).*(1+x(:,1)).*cos(pi*x(:,2))];
        data.boundtype=[-1 -1 -2 -2];         % -1 inflow, -2 outflow
        data.Coordinates=[-0.125 -0.125; 1 -0.125; 1 1; -0.125 1];
        data.Elements=[1 2 4;2 3 4];
        data.BiForms=[-1.3069096523286656, 16.49280590717239, 0, -0.7572149944481313;...
            2.6138193184362004, -1.3069096523286656, 0, -0.17649400565242843];
            
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
        % - curl u
        data.CURL_U_EX_Handle=@(x,varargin)(pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1))-pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2)));
        
        data.boundtype=[-1 -1 -2 -2];
        data.Coordinates=[-1 -1; 1 -1; 1 1; -1 1];
        data.Elements=[1 2 4;2 3 4];
        data.BiForms=[-4*pi^2, 8*pi^2, 0, 0 ; ...
                               4*pi^2, -4*pi^2, 0 ,0];
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
        % - curl u
        data.CURL_U_EX_Handle=@(x,varargin)zeros(size(x,1),1);            
        data.boundtype=[-1 -1 -2 -2];
        data.Coordinates=[-1 -1; 1 -1; 1 1; -1 1];
        data.Elements=[1 2 4;2 3 4];
        data.BiForms=[0, 32/3,0,-88/15; 0, 0, 0, 0];
    
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
        % - curl u
        data.CURL_U_EX_Handle=@(x,varargin)zeros(size(x,1),1) ;
        data.boundtype=[-1 -1 -2 -2];
        data.Coordinates=[-1 -1; 1 -1; 1 1; -1 1];
        data.Elements=[1 2 4;2 3 4];
        data.BiForms=[0,8,0,0; 0,0,0,0];
    
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
                                       (2+sin(pi.*x(:,1))).*(-pi.*(-1+x(:,1).^2).*cos(pi.*x(:,2))-2.*x(:,2).*sin(pi.*x(:,1)))];
        % curl curl u
        data.SOL3_Handle=@(x,varargin)[2.*(-pi.*x(:,1).*cos(pi.*x(:,2))+sin(pi.*x(:,1))) ...
                                       2.*(-pi.*x(:,2).*cos(pi.*x(:,1))+sin(pi.*x(:,2)))];
        % - curl u
        data.CURL_U_EX_Handle=@(x,varargin)-2.*x(:,2).*sin(pi.*x(:,1))+2.*x(:,1).*sin(pi.*x(:,2));
        
        %data.Coordinates=[-1 -1; 1 -1; 1 1; -1 1];  
        data.Coordinates=[-0.125 -0.125; 1 -0.125; 1 1; -0.125 1];
        %data.Elements=[1 2 4;2 3 4];
        %data.boundtype=[-1 -1 -2 -2];
        data.Elements=[1 2 3;1 3 4];
        data.boundtype=[-1 -1 -2 -2];
   
    case 6
        % velocity v
        data.V_Handle=@(x,varargin)[-ones(size(x,1),1)-x(:,2) ones(size(x,1),1)+x(:,1)];
        % solution u
        data.U_EX_Handle=@(x,varargin)[(1-x(:,2)).*(1+x(:,2)).*sin(pi*x(:,1)) ...
                                       (1-x(:,1)).*(1+x(:,1)).*sin(pi*x(:,2))];
        %  v *u                  
        data.V_U_EX_Handle=@(x,varargin)(-1-x(:,2)).*(1-x(:,2).^2).*sin(pi*x(:,1))+...
                                        (1+x(:,1)).*(1-x(:,1).^2).*sin(pi*x(:,2)) ;
        %  v x curl u                          
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
        data.BiForms=[-8-96/pi^4+64/pi^2, 64/105*(14-945/pi^4+4*pi^2), 0, -(16/pi^3) ; ...
            704/45+384/pi^4-392/(3*pi^2), -8-96/pi^4+64/pi^2, 0, -((16*(-15+pi^2))/(15*pi^3))];
        
    case 7
        % velocity v
        data.V_Handle=@(x,varargin)[ones(size(x,1),1) ones(size(x,1),1)];
        % solution u
        data.U_EX_Handle=@(x,varargin)[1-x(:,2) 1+x(:,1)];
        %  v *u                  
        data.V_U_EX_Handle=@(x,varargin)2+x(:,1)-x(:,2);
        %  - v x curl u                          
        data.SOL1_Handle=@(x,varargin)[-2*ones(size(x,1),1) 2*ones(size(x,1),1)];
        % grad(v*u)
        data.SOL2_Handle=@(x,varargin)[ones(size(x,1),1) -ones(size(x,1),1)];
        % curl curl u
        data.SOL3_Handle=@(x,varargin)zeros(size(x,1),2)  ;
        % -curl u
        data.CURL_U_EX_Handle=@(x,varargin)-2*ones(size(x,1),1) ;
        data.boundtype=[-1 -1 -2 -2];
        data.Coordinates=[-1 -1; 1 -1; 1 1; -1 1];
        data.Elements=[1 2 4;2 3 4];   
        data.BiForms=[-16,8,0,0; 32,-16,0,0];
        
    case 8
        % velocity v
        data.V_Handle=@(x,varargin)ones(size(x,1),2);
        % solution u
        data.U_EX_Handle=@(x,varargin)[-sin(pi.*x(:,2)).*sin(pi.*x(:,1)) ...
                                       sin(pi.*x(:,1)).*sin(pi.*x(:,2))];
        %  v *u                  
        data.V_U_EX_Handle=@(x,varargin)zeros(size(x,1),1);
        %  - v x curl u                          
        data.SOL1_Handle=@(x,varargin)[-pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1))-pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2)) ...
                                                        +pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1))+pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2))]; 
        % grad(v*u)
        data.SOL2_Handle=@(x,varargin)zeros(size(x,1),2); 
        % curl curl u
        data.SOL3_Handle=@(x,varargin)[pi.^2.*cos(pi.*x(:,2)).*cos(pi.*x(:,1))-pi.^2.*sin(pi.*x(:,1)).*sin(pi.*x(:,2)) ...
                                      -pi.^2.*cos(pi.*x(:,2)).*cos(pi.*x(:,1))+pi.^2.*sin(pi.*x(:,1)).*sin(pi.*x(:,2))];
        % - curl u
        data.CURL_U_EX_Handle=@(x,varargin)(-pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1))-pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2)));
        
        data.boundtype=[-1 -1 -2 -2];
        data.Coordinates=[-1 -1; 1 -1; 1 1; -1 1];
        data.Elements=[1 2 4;2 3 4];
        data.BiForms=[-4*pi^2, 8*pi^2, 0, 0 ; ...
                               4*pi^2, -4*pi^2, 0 ,0];
                           
       case 9
        % velocity v
        data.V_Handle=@(x,varargin)[ones(size(x,1),1) 0.5.*ones(size(x,1),1)];
        % solution u
        data.U_EX_Handle=@(x,varargin)...
                       - [4/5*(x(:,2).^2-1).*sin(pi*x(:,1)) + 2/5*(x(:,1).^2-1).*sin(pi*x(:,2)) ...
                         2/5*(x(:,2).^2-1).*sin(pi*x(:,1)) + 1/5*(x(:,1).^2-1).*sin(pi*x(:,2))];
        %  v *u                  
        data.V_U_EX_Handle=@(x,varargin)...
                       - ((x(:,2).^2-1).*sin(pi*x(:,1)) + 1/2*(x(:,1).^2-1).*sin(pi*x(:,2)));
        %  - v x curl u                          
        data.SOL1_Handle=@(x,varargin)...
                       1/5*[ pi*(x(:,2).^2-1).*cos(pi*x(:,1)) - pi*(x(:,1).^2-1).*cos(pi*x(:,2))-...
                              4*x(:,2).*sin(pi*x(:,1)) + x(:,1).*sin(pi*x(:,2)) ...
                            -2*(pi*(x(:,2).^2-1).*cos(pi*x(:,1)) - pi*(x(:,1).^2-1).*cos(pi*x(:,2))-...
                              4*x(:,2).*sin(pi*x(:,1)) + x(:,1).*sin(pi*x(:,2)))];    % - v x curl u
        % grad(v*u)
        data.SOL2_Handle=@(x,varargin)...
                        [-pi*(x(:,2).^2-1).*cos(pi*x(:,1))-x(:,1).*sin(pi*x(:,2)) ...
                         -pi/2*(x(:,1).^2-1).*cos(pi*x(:,2))-2*x(:,2).*sin(pi*x(:,1))]; % grad (v u)
        % curl curl u
        data.SOL3_Handle=@(x,varargin)...                  
                  1/5*[8*sin(pi.*x(:,1))-2*pi*( 2*x(:,2).*cos(pi*x(:,1))+x(:,1).*cos(pi*x(:,2))+pi*(x(:,1).^2-1).*sin(pi*x(:,2))) ...
                          2*sin(pi.*x(:,2))-2*pi*( 4*x(:,2).*cos(pi*x(:,1))+2*x(:,1).*cos(pi*x(:,2))+pi*(x(:,2).^2-1).*sin(pi*x(:,1)))]; % curl curl u
        % - curl u
        data.CURL_U_EX_Handle=@(x,varargin) ...
                         -(-2/5*pi*(x(:,2).^2-1).*cos(pi*x(:,1)) + 2/5*pi*(x(:,1).^2-1).*cos(pi*x(:,2))+...
                            8/5*x(:,2).*sin(pi*x(:,1)) -2/5* x(:,1).*sin(pi*x(:,2)));    % curl u
        
        % div u
        data.Div_U_EX_Handle=@(x,varargin)[pi*(1-x(:,2)).*(1+x(:,2)).*cos(pi*x(:,1))+...
                                pi*(1-x(:,1)).*(1+x(:,1)).*cos(pi*x(:,2))];
        data.boundtype=[-1 -1 -2 -2];         % -1 inflow, -2 outflow
        data.Coordinates=[-1 -1; 1 -1; 1 1; -1 1];
        data.Elements=[1 2 4;2 3 4];
        data.BiForms=[-1.3069096523286656, 16.49280590717239, 0, -0.7572149944481313;...
            2.6138193184362004, -1.3069096523286656, 0, -0.17649400565242843];
       
     case 10
        % velocity v 
        data.V_Handle=@(x,varargin)[ones(size(x,1),1) 0.5.*ones(size(x,1),1)];
        % solution u
        data.U_EX_Handle=@(x,varargin)...
                       - [4/5*(x(:,2).^2-1).*sin(pi*x(:,1)) + 2/5*(x(:,1).^2-1).*sin(pi*x(:,2)) ...
                         2/5*(x(:,2).^2-1).*sin(pi*x(:,1)) + 1/5*(x(:,1).^2-1).*sin(pi*x(:,2))]-...
                         [(1-x(:,2)).*(1+x(:,2)).*sin(pi*x(:,1)) ...
                          (1-x(:,1)).*(1+x(:,1)).*sin(pi*x(:,2))];
        %  v *u                  
        data.V_U_EX_Handle=@(x,varargin)(1-x(:,2)).*(1+x(:,2)).*sin(pi*x(:,1))+...
                                                    0.5*(1-x(:,1)).*(1+x(:,1)).*sin(pi*x(:,2))+...
                                                    ((x(:,2).^2-1).*sin(pi*x(:,1)) + 1/2*(x(:,1).^2-1).*sin(pi*x(:,2)));
        %  - v x curl u                          
        data.SOL1_Handle=@(x,varargin)...
                       1/5*[ pi*(x(:,2).^2-1).*cos(pi*x(:,1)) - pi*(x(:,1).^2-1).*cos(pi*x(:,2))-...
                              4*x(:,2).*sin(pi*x(:,1)) + x(:,1).*sin(pi*x(:,2)) ...
                            -2*(pi*(x(:,2).^2-1).*cos(pi*x(:,1)) - pi*(x(:,1).^2-1).*cos(pi*x(:,2))-...
                              4*x(:,2).*sin(pi*x(:,1)) + x(:,1).*sin(pi*x(:,2)))]-...
                              [-x(:,2).*sin(pi.*x(:,1))+x(:,1).*sin(pi.*x(:,2)) ...
                               2.*x(:,2).*sin(pi.*x(:,1))-2.*x(:,1).*sin(pi.*x(:,2))];    % - v x curl u
        % grad(v*u)
        data.SOL2_Handle=@(x,varargin)...
                        [-pi*(x(:,2).^2-1).*cos(pi*x(:,1))-x(:,1).*sin(pi*x(:,2)) ...
                         -pi/2*(x(:,1).^2-1).*cos(pi*x(:,2))-2*x(:,2).*sin(pi*x(:,1))]-...
                        [-pi*(x(:,2).^2-1).*cos(pi*x(:,1))-x(:,1).*sin(pi*x(:,2)) ...
                         -pi/2*(x(:,1).^2-1).*cos(pi*x(:,2))-2*x(:,2).*sin(pi*x(:,1))];% grad (v u)
        % curl curl u
        data.SOL3_Handle=@(x,varargin)...                  
                  1/5*[8*sin(pi.*x(:,1))-2*pi*( 2*x(:,2).*cos(pi*x(:,1))+x(:,1).*cos(pi*x(:,2))+pi*(x(:,1).^2-1).*sin(pi*x(:,2))) ...
                          2*sin(pi.*x(:,2))-2*pi*( 4*x(:,2).*cos(pi*x(:,1))+2*x(:,1).*cos(pi*x(:,2))+pi*(x(:,2).^2-1).*sin(pi*x(:,1)))]-...
                          [2.*(-pi.*x(:,1).*cos(pi.*x(:,2))+sin(pi.*x(:,1))) ...
                           2.*(-pi.*x(:,2).*cos(pi.*x(:,1))+sin(pi.*x(:,2)))]; % curl curl u
        % - curl u
        data.CURL_U_EX_Handle=@(x,varargin) ...
                         -(-2/5*pi*(x(:,2).^2-1).*cos(pi*x(:,1)) + 2/5*pi*(x(:,1).^2-1).*cos(pi*x(:,2))+...
                            8/5*x(:,2).*sin(pi*x(:,1)) -2/5* x(:,1).*sin(pi*x(:,2)))-...
                          (-2.*x(:,2).*sin(pi.*x(:,1))+2.*x(:,1).*sin(pi.*x(:,2)));    % curl u
        
        % div u
        data.Div_U_EX_Handle=@(x,varargin)[pi*(1-x(:,2)).*(1+x(:,2)).*cos(pi*x(:,1))+...
                                pi*(1-x(:,1)).*(1+x(:,1)).*cos(pi*x(:,2))];
        data.boundtype=[-1 -1 -2 -2];         % -1 inflow, -2 outflow
        data.Coordinates=[-1 -1; 1 -1; 1 1; -1 1];
        data.Elements=[1 2 4;2 3 4];
        data.BiForms=[-1.3069096523286656, 16.49280590717239, 0, -0.7572149944481313;...
            2.6138193184362004, -1.3069096523286656, 0, -0.17649400565242843];
        
        case 11
        % velocity v
        data.V_Handle=@(x,varargin)[ones(size(x,1),1) 0.5.*ones(size(x,1),1)];
        
        %  solution u
        data.U_EX_Handle=@(x,varargin)[sin(pi.*x(:,2)).*sin(pi.*x(:,1)) ...
                                                       sin(pi.*x(:,1)).*sin(pi.*x(:,2))];
        %  v*u                  
        data.V_U_EX_Handle=@(x,varargin)(sin(pi.*x(:,2)).*sin(pi.*x(:,1))+0.5*sin(pi.*x(:,1)).*sin(pi.*x(:,2)));
        
        %  - v x curl u                          
        data.SOL1_Handle=@(x,varargin)[+0.5*pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1))-0.5*pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2)) ...
                                                        -pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1))+pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2))]; 
        %  grad(v*u)
        data.SOL2_Handle=@(x,varargin)[1.5.*pi.*sin(pi.*x(:,2)).*cos(pi.*x(:,1))...
                                                       1.5.*pi.*sin(pi.*x(:,1)).*cos(pi.*x(:,2))]; 
        %  curl curl u
        data.SOL3_Handle=@(x,varargin)[pi.^2.*cos(pi.*x(:,2)).*cos(pi.*x(:,1))+pi.^2.*sin(pi.*x(:,1)).*sin(pi.*x(:,2)) ...
                                                       +pi.^2.*cos(pi.*x(:,2)).*cos(pi.*x(:,1))+pi.^2.*sin(pi.*x(:,1)).*sin(pi.*x(:,2))];
        %  - curl u
        data.CURL_U_EX_Handle=@(x,varargin)(pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1))-pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2)));
        
        data.boundtype=[-1 -1 -2 -2];
        data.Coordinates=[-1.5 -1.5; 1.5 -1.5; 1.5 1.5; -1.5 1.5];
        data.Elements=[1 2 4;2 3 4];
        data.BiForms=[-4*pi^2, 8*pi^2, 0, 0 ; ...
                               4*pi^2, -4*pi^2, 0 ,0];
                           
    case 12
        % velocity v
        l1=0.5; 
        l2=0.5;
        data.V_Handle=@(x,varargin)[l1*ones(size(x,1),1) l2*ones(size(x,1),1)];
        
        %  solution u
        data.U_EX_Handle=@(x,varargin)[sin(pi.*x(:,2)).*sin(pi.*x(:,1)) ...
                                                      sin(pi.*x(:,1)).*sin(pi.*x(:,2))];
        %  v*u                  
        data.V_U_EX_Handle=@(x,varargin)(l1*sin(pi.*x(:,2)).*sin(pi.*x(:,1))+l2*sin(pi.*x(:,1)).*sin(pi.*x(:,2)));
        
        %  - v x curl u                          
        data.SOL1_Handle=@(x,varargin)[+l2*pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1))-l2*pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2)) ...
                                                       -l1*pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1))+l1*pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2))]; 
        %  grad(v*u)
        data.SOL2_Handle=@(x,varargin)[(l1+l2).*pi.*sin(pi.*x(:,2)).*cos(pi.*x(:,1))...
                                                       (l1+l2).*pi.*sin(pi.*x(:,1)).*cos(pi.*x(:,2))]; 
        %  curl curl u
        data.SOL3_Handle=@(x,varargin)[pi.^2.*cos(pi.*x(:,2)).*cos(pi.*x(:,1))+pi.^2.*sin(pi.*x(:,1)).*sin(pi.*x(:,2)) ...
                                                      +pi.^2.*cos(pi.*x(:,2)).*cos(pi.*x(:,1))+pi.^2.*sin(pi.*x(:,1)).*sin(pi.*x(:,2))];
        %  - curl u
        data.CURL_U_EX_Handle=@(x,varargin)(pi.*cos(pi.*x(:,2)).*sin(pi.*x(:,1))-pi.*cos(pi.*x(:,1)).*sin(pi.*x(:,2)));
        
        data.boundtype=[-1 -1 -2 -2];
        data.Coordinates=[-1.5 -1.5; 1.5 -1.5; 1.5 1.5; -1.5 1.5];
        data.Elements=[1 2 3;1 3 4];
        
    case 13
        % velocity v
        data.V_Handle=@(x,varargin)[ones(size(x,1),1) ones(size(x,1),1)];
        % solution u
        data.U_EX_Handle=@(x,varargin)[1+x(:,1) 1+x(:,2)];
        %  v *u                  
        data.V_U_EX_Handle=@(x,varargin)2+x(:,1)+x(:,2);
        %  - v x curl u                          
        data.SOL1_Handle=@(x,varargin)zeros(size(x,1),2);
        % grad(v*u)
        data.SOL2_Handle=@(x,varargin)[ones(size(x,1),1) ones(size(x,1),1)];
        % curl curl u
        data.SOL3_Handle=@(x,varargin)zeros(size(x,1),2);
        % -curl u
        data.CURL_U_EX_Handle=@(x,varargin)zeros(size(x,1),1);
        data.boundtype=[-1 -1 -2 -2];
        data.Coordinates=[-1 -1; 1 -1; 1 1; -1 1];
        data.Elements=[1 2 4;2 3 4];   
        data.BiForms=[-16,8,0,0; 32,-16,0,0];
        
      case 14
        % velocity v
        data.V_Handle=@(x,varargin)[zeros(size(x,1),1) ones(size(x,1),1)];
        % solution u
        data.U_EX_Handle=@(x,varargin)zeros(size(x,1),2);
        %  v *u                  
        data.V_U_EX_Handle=@(x,varargin)zeros(size(x,1),1);
        %  - v x curl u                          
        data.SOL1_Handle=@(x,varargin)zeros(size(x,1),2);
        % grad(v*u)
        data.SOL2_Handle=@(x,varargin)zeros(size(x,1),2);
        % curl curl u
        data.SOL3_Handle=@(x,varargin)zeros(size(x,1),2);
        % -curl u
        data.CURL_U_EX_Handle=@(x,varargin)zeros(size(x,1),1);
        data.boundtype=[-1 -3 -4 -2];
        data.Coordinates=[-1 -1; 1 -1; 1 1; -1 1];
        data.Elements=[1 2 4;2 3 4];   
        data.BiForms=[-16,8,0,0; 32,-16,0,0];
        
    case 15
        H1=@(x) sin(pi.*x(:,1)).*(1-x(:,2).^2);
        D1H1=@(x) pi*(1-x(:,2).^2).*cos(pi.*x(:,1));
        D2H1=@(x) -2*x(:,2).*sin(pi.*x(:,1));
        D11H1=@(x) pi^2*(1-x(:,2).^2).*sin(pi.*x(:,1));
        D12H1=@(x) -2*pi*x(:,2).*cos(pi.*x(:,1));
        D22H1=@(x) -2*sin(pi.*x(:,1));

        H2=@(x) (1-x(:,2).^2).*(1-x(:,1).^2);
        D1H2=@(x)(1-x(:,2).^2).*(-2*x(:,1));
        D2H2=@(x)(-2*x(:,2)).*(1-x(:,1).^2);
        D11H2=@(x)(1-x(:,2).^2).*(-2);
        D12H2=@(x)(-2*x(:,2)).*(-2*x(:,1));
        D22H2=@(x)(-2).*(1-x(:,1).^2);

        H1=@(x)(-2*x(:,2)).*(1-x(:,1).^2);
        D1H1=@(x)(-2*x(:,2)).*(-2*x(:,1));
        D2H1=@(x)-2*(1-x(:,1).^2);
        D11H1=@(x)(-2*x(:,2)).*(-2);
        D12H1=@(x)(-2)*(-2*x(:,1));
        D22H1=@(x)zeros(size(x,1),1);
        H2=@(x) (1-x(:,2).^2).*(2*x(:,1));
        D1H2=@(x)(1-x(:,2).^2)*2;
        D2H2=@(x)(-2*x(:,2)).*(2*x(:,1));
        D11H2=@(x)zeros(size(x,1),1);
        D12H2=@(x)(-2*x(:,2))*2;
        D22H2=@(x)(-2).*(2*x(:,1));

        %        v1=0.5;
        %        v2=0.5;
        %         H1 = @(x) 1 * (1-x(:,1).^2).^2 .*(x(:,2)-x(:,2).^3);
        %         D1H1=@(x) -4 * (x(:,1)-x(:,1).^3) .*(x(:,2)-x(:,2).^3);
        %         D2H1=@(x) 1 * (1-x(:,1).^2).^2 .*(1-3*x(:,2).^2);
        %         D11H1=@(x) -4 * (1-3*x(:,1).^2) .*(x(:,2)-x(:,2).^3);
        %         D12H1=@(x) -4 * (x(:,1)-x(:,1).^3) .*(1-3*x(:,2).^2);
        %         D22H1=@(x) -6 * (1-x(:,1).^2).^2 .*x(:,2);
        %         H2=@(x) -1 * (x(:,1)-x(:,1).^3).*(1-x(:,2).^2).^2 ;
        %         D1H2=@(x) -1 * (1-3*x(:,1).^2).*(1-x(:,2).^2).^2;
        %         D2H2=@(x) 4 * (x(:,1)-x(:,1).^3).*(x(:,2)-x(:,2).^3);
        %         D11H2=@(x) 6 * x(:,1).*(1-x(:,2).^2).^2;
        %         D12H2=@(x) 4 * (1-3*x(:,1).^2).*(x(:,2)-x(:,2).^3);
        %         D22H2=@(x) 4 * (x(:,1)-x(:,1).^3).*(1-3*x(:,2).^2);
        %
        % %         v1=0.5;
        % %         v2=0.5;

        DIR1 = @(x)  (1-x(:,1).^2).^2 .*(x(:,2)-x(:,2).^3);
        D1DIR1=@(x) -4 * (x(:,1)-x(:,1).^3) .*(x(:,2)-x(:,2).^3);
        D2DIR1=@(x)  (1-x(:,1).^2).^2 .*(1-3*x(:,2).^2);
        DIR2=@(x) -1 * (x(:,1)-x(:,1).^3).*(1-x(:,2).^2).^2 ;
        D1DIR2=@(x) -1 * (1-3*x(:,1).^2).*(1-x(:,2).^2).^2;
        D2DIR2=@(x) 4 * (x(:,1)-x(:,1).^3).*(x(:,2)-x(:,2).^3);

        %         v1=0.66;
        %         v2=1;
        %         DIR1=@(x)v1*(1-x(:,2).^2).*(1-x(:,1).^2);
        %         %DIR1=@(x)v1*ones(size(x,1),1);
        %         D1DIR1=@(x)v1*(1-x(:,2).^2).*(-2*x(:,1));
        %         D2DIR1=@(x)v1*(-2*x(:,2)).*(1-x(:,1).^2);
        %         DIR2=@(x)sin(pi.*x(:,2)).*sin(pi.*x(:,1));
        %         %DIR2=@(x)v2*ones(size(x,1),1);
        %         D1DIR2=@(x)pi*sin(pi.*x(:,2)).*cos(pi.*x(:,1));
        %         D2DIR2=@(x)pi*cos(pi.*x(:,2)).*sin(pi.*x(:,1));

        DIR1 = @(x)  1*ones(size(x(:,1)));
        D1DIR1 = @(x)  0*ones(size(x(:,1)));
        D2DIR1 = @(x)  0*ones(size(x(:,1)));
        DIR2=@(x) 1*ones(size(x(:,1)));
        D1DIR2=@(x) 0*ones(size(x(:,1)));
        D2DIR2=@(x) 0*ones(size(x(:,1)));

        data.V_Handle=@(x,varargin)[DIR1(x), DIR2(x)];
        data.JacV_Handle=@(x,varargin)[D1DIR1(x) D1DIR2(x) D2DIR1(x) D2DIR2(x)];
        data.U_EX_Handle=@(x,varargin)[H1(x) H2(x)];

        %  - v x curl u                          
        data.SOL1_Handle=@(x,varargin)...
            [-DIR2(x).*(D1H2(x)-D2H1(x)) ...
            DIR1(x).*(D1H2(x)-D2H1(x))];

        %  v *u                 
        data.V_U_EX_Handle=@(x,varargin)...
            [DIR1(x).*H1(x)+DIR2(x).*H2(x)];

        % grad(v*u)
        data.SOL2_Handle=@(x,varargin)...
            [D1DIR1(x).*H1(x)+DIR1(x).*D1H1(x)+D1DIR2(x).*H2(x)+DIR2(x).*D1H2(x) ...
             D2DIR1(x).*H1(x)+DIR1(x).*D2H1(x)+D2DIR2(x).*H2(x)+DIR2(x).*D2H2(x)];
        
        % curl curl u
        data.SOL3_Handle=@(x,varargin)...
            [-D22H1(x) + D12H2(x) ...
               D12H1(x) - D11H2(x)];
        % -curl u
        data.CURL_U_EX_Handle=@(x,varargin)-D1H2(x)+D2H1(x);
        %        data.boundtype=[-1 -1 -1 -1];
        data.Coordinates=[-1 -1; 1 -1; 1 1; -1 1];
        %data.Coordinates=[0 0; 1 0; 1 1; 0 1];
        %data.Coordinates=[0.25 0.25 ; 0.6 0.25; 0.6 0.6; 0.25 0.6];
        data.Elements=[1 2 4;2 3 4];   
        data.boundtype=[-1 -1 -2 -2];     
        %data.Elements=[1 2 3;1 3 4];   
        %data.boundtype=[-1 -1 -2 -2];     
        
     end
                 
d=data;
return