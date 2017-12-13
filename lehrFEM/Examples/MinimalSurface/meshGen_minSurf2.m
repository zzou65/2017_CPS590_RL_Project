function meshGen_minSurf2
% meshGen_minSurf2
%   Run script for creating graded square meshes for the minimal surface
%   problem.

% Copyright 2006-2006 Kari Borset
% SAM - Seminar for Applied Mathematics
% ETH-Zentrum
% CH-8092 Zurich, Switzerland

clc

% Initialize constants

h_min = 0.001;            % Smallest mesh width
% A = load('gradFunc.mat');
% g = A.o(:,2);               % grading function
% h_grad = 1/(1+size(g,1));   % stepzise used in gradingfunction
g_func = @(x,varargin)abs(x).^2.1;
g = g_func(linspace(0,1,10000))';
DF = 15;                    % Number of gridpoints in x-direction (in one half plane) for first mesh
N = 5;                      % Number of meshes
% HHANDLE = @(x,varargin)1+4.3842*x.^5-3.1918*x.^4-5.18*x.^3+3.4998*x.^2-0.68*x;   % Grading in x-direction
%HHANDLE = @(x,varargin)sign(x)*g(min(floor(abs(x)/h_grad)+1,size(g)));   % Grading in x-direction


for i = 1:N,
    t = linspace(0,1,DF/2);    
    x = HHANDLE(t,g);
    
    h_max = norm(x(1)-x(2));
    for j=1:5,
        h_max_b = max(norm(x(i)-x(i+1)),h_max);
    end
    for j=1:4,
        if 2*x(1)<h_max,
            x = x(2:end);
            t = t(2:end);
        end
    end
    x = [-x x(1:end-1)];
    t = [-t t(2:end)];
   % plot(t,x,'rx')
%     figure()
%     plot(t,zeros(size(t)),'ro', ...
%          x,ones(size(x)),'bo')
     
        
%     figure()
%     plot(x,zeros(size(x)),'ro')
    
    y = linspace(-1,1,2*DF);

    [X Y] = meshgrid(x,y);
    Mesh.Coordinates = [X(:) Y(:)];
    Mesh = add_Elements(Mesh);
    save_Mesh(Mesh, ...
              ['Coord_Sqr_' int2str(i) '.dat'], ...
              ['Elem_Sqr_' int2str(i) '.dat']);
    DF = floor(DF*1.5);
%     plot_Qual(Mesh)
end


function f = HHANDLE(x,g,varargin)

n = length(x);
h = 1/(size(g,1));

s = sign(x);
x = abs(x);
for i = 1:n,
    if x(i)==0,
        f(i)=1;
    elseif x(i)==1,
        f(i)=0;
    else
        inda = floor(x(i)/h)+1; a = inda*h;
        indb = floor(x(i)/h)+2; b = indb*h;
        f(i) = 1-((b-x(i))*g(inda) + (x(i)-a)*g(indb))/(b-a);
    end
end
f = f.*s;

return

% function f = HHANDLE(x,varargin)
% if x<0.35,
%     s = 1;
% elseif x>0.8,
%     s = 3;
% else
%     s = 2;
% end
% switch s
%     case 1
%         f = 96.63*x.^4-332.76*x.^3+428.59*x.^2-246.48*x+54.03;
%     case 2
%         f = -209.64*x.^6+486.67*x.^5-463.51*x.^4+229.81*x.^3-62.82*x.^2+8.93*x+0.49;
%     case 3
%         f =  1-0.072*x.^2+0.01;
% end
% if x == 0,
%     f=0;
% elseif x==1,
%     f=1;
% end
% return
% 
% 
% for i = 1:N,
%     Mesh = Sqr_graded(h);
%     save_Mesh(Mesh, ...
%               ['Coord_Sqr_' int2str(i) '.dat'], ...
%               ['Elem_Sqr_' int2str(i) '.dat']);
%     close all
%     h = h/2;
% end