% plots the vector valued linear shape functions on the reference triangle

%   Copyright 2009-2009 Christoph Wiesmeyr
%   SAM - Seminar for Applied Mathematics
%   ETH-Zentrum
%   CH-8092 Zurich, Switzerland

x = linspace(0,1,10);
y = linspace(0,1,10);



count = 0;

for i = 1:length(x)
    for j = 1:length(y)
        if x(i) + y(j) <= 1
            count = count+1;
            meshpx(count) = x(i);
            meshpy(count) = y(j);
            Shap = shap_LFE2([x(i) y(j)]);
            u1(count) = Shap(1);
            v1(count) = Shap(2);
            u2(count) = Shap(3);
            v2(count) = Shap(4);
            u3(count) = Shap(5);
            v3(count) = Shap(6);
            u4(count) = Shap(7);
            v4(count) = Shap(8);
            u5(count) = Shap(9);
            v5(count) = Shap(10);
            u6(count) = Shap(11);
            v6(count) = Shap(12);
        end
    end
end


% first shape function
figure
hold on
quiver(meshpx,meshpy,u1,v1);
plot([0 0 0; 0 1 1], [0 0 1; 1 0 0],'r');
hold off


% second shape function
figure
hold on
quiver(meshpx,meshpy,u2,v2);
plot([0 0 0; 0 1 1], [0 0 1; 1 0 0],'r');
hold off

% third shape function
figure
hold on
quiver(meshpx,meshpy,u3,v3);
plot([0 0 0; 0 1 1], [0 0 1; 1 0 0],'r');
hold off

% 4th shape function
figure
hold on
quiver(meshpx,meshpy,u4,v4);
plot([0 0 0; 0 1 1], [0 0 1; 1 0 0],'r');
hold off


% 5th shape function
figure
hold on
quiver(meshpx,meshpy,u5,v5);
plot([0 0 0; 0 1 1], [0 0 1; 1 0 0],'r');
hold off

% 6th shape function
figure
hold on
quiver(meshpx,meshpy,u6,v6);
plot([0 0 0; 0 1 1], [0 0 1; 1 0 0],'r');
hold off


clear all;