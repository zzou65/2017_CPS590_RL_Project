function d=getData(n)

switch n
    %
    case 1 % easy smooth
        % velocity v
        data.V_Handle=...
            @(x,varargin)ones(size(x,1),2);
        % solution u
        data.U_EX_Handle=...
            @(x,varargin)sin(pi.*x(:,1))+cos(pi.*x(:,2));
        % solution grad u
        data.GRAD_U_EX_Handle=...
            @(x,varargin)[pi.*cos(pi.*x(:,1)) -pi.*sin(pi.*x(:,2))];                    
        % -a*div grad u+b*v grad u
        data.SOL_Handle=...
            @(x,elemflag,a)a*(pi.^2.*sin(pi.*x(:,1))+pi.^2.*cos(pi.*x(:,2)))+...
            (pi.*cos(pi.*x(:,1))-pi.*sin(pi.*x(:,2)));
        data.Coordinates=[-1 -1; 1 -1; 1 1; -1 1];
        data.Elements=[1 2 4;2 3 4];
        data.boundtype=[-1 -1 -1 -1];
        
     case 2 % triangle domain, boundary layer
        % velocity v
        data.V_Handle=...
            @(x,varargin)[ones(size(x,1),1) zeros(size(x,1),1)];
        % solution u
        data.U_EX_Handle=...
            @(x,d,a)x(:,1)-1./(1-exp(-1/a)).*(exp((x(:,1)-1)/a)-exp(-1/a));
        % solution grad u
        data.GRAD_U_EX_Handle=...
            @(x,varargin)0;                    
        % -a*div grad u+b*v grad u
        data.SOL_Handle=@(x,elemflag,a,b)ones(size(x,1),1);
        data.Coordinates=[0 0; 1 -1; 1  1];   
        data.Elements= [1 2 3];
        data.boundtype=[-1 -1 -1];    
        
      case 3 % layer at diagonal
        % velocity v
        data.V_Handle=...
            @(x,varargin)[ones(size(x,1),1) ones(size(x,1),1)];
        % solution u
        data.U_EX_Handle=@(x,d,a)(x(:,2)<=x(:,1));
        % solution grad u
        data.GRAD_U_EX_Handle=...
            @(x,varargin)0;                    
        % -a*div grad u+b*v grad u
        data.SOL_Handle=@(x,elemflag,a,b)zeros(size(x,1),1);
        data.Coordinates=[0 0; 1 0; 1 1; 0 1];
        data.Elements=[1 2 3;1 3 4];
        data.boundtype=[-1 -1 -1 -1];   
      
     case 4  % Houston, Sueli, Schwab hp-FEM... 2002
        % velocity v
        data.V_Handle=...
            @(x,varargin)[ones(size(x,1),1) ones(size(x,1),1)];
        % solution u
        data.U_EX_Handle=@(x,d,a)...
            x(:,1)+x(:,2).*(1-x(:,1))+(exp(-1/a)-exp(-(1-x(:,1)).*(1-x(:,2))/a))/(1-exp(-1/a));
        % solution grad u
        data.GRAD_U_EX_Handle=...
            @(x,varargin)0;                    
        % -a*div grad u+b*v grad u
        data.SOL_Handle=@(x,elemflag,a,b) ...
            (2-x(:,1)-x(:,2))+(exp((x(:,1)+x(:,2)-x(:,1).*x(:,2))/a).*...
            ((-1+x(:,1)).*x(:,1)+(-1+x(:,2)).*x(:,2)))/(a*(-1+exp(1/a)));
        
        data.Coordinates=[0 0; 1 0; 1 1; 0 1];
        data.Elements=[1 2 3;1 3 4];
        data.boundtype=[-1 -1 -1 -1];     
        
      case 5 % Mathiess Dipl
        % velocity v
        data.V_Handle=...
            @(x,varargin)[2*ones(size(x,1),1) 3*ones(size(x,1),1)];
        % solution u
        data.U_EX_Handle=@(x,d,a)...
            x(:,1).*x(:,2).^2-x(:,2).^2.*exp(2*(x(:,1)-1)/a)-x(:,1).*exp(3*(x(:,2)-1)/a)...
            +exp((2*(x(:,1)-1)/a+3*(x(:,2)-1))/a);
        % solution grad u
        data.GRAD_U_EX_Handle=...
            @(x,varargin)0;                    
        % -a*div grad u+b*v grad u
        data.SOL_Handle=@(x,elemflag,a,b)...
            2*(-exp(3*(-1+x(:,2))/a)+a*(exp(2*(-1+x(:,1))/a)-x(:,1))...
            -3*exp(2/a*(-1+x(:,1))).*x(:,2)+x(:,2).*(3*x(:,1)+x(:,2)));
        data.Coordinates=[0 0; 1 0; 1 1; 0 1];
        data.Elements=[1 2 3;1 3 4];
        data.boundtype=[-1 -1 -1 -1];     
      
     case 6 %  Franca  CMAME 92
        % velocity v
        data.V_Handle=...
            @(x,varargin)[-x(:,2) x(:,1)];
        % solution u
        data.U_EX_Handle=@(x,d,a)...
            1/2*(cos(4*pi*sqrt(x(:,1).^2+x(:,2).^2)+pi)+1);
        % solution grad u
        data.GRAD_U_EX_Handle=...
            @(x,varargin)0;                    
        % -a*div grad u+b*v grad u
        data.SOL_Handle=@(x,elemflag,a,b)0;
        data.Coordinates=[0 0; 0 -0.5 ; 0.5 -0.5; ...
                                   0.5 0; 0.5 0.5; 0 0.5; ...
                                 -0.5 0.5; -0.5 0; -0.5 -0.5]; 
        data.Elements=[1 2 3; 3 4 1; 5 6 4; 1 4 6; ...
                               1 6 7; 1 7 8;];
        data.boundtype=[-1 -2 -1 -1 -1 -1 -1 -1];       
        
end
        
        
d=data;
return