clc
clear
load("examples\mouse.mat");
u = 0:0.01:1;
figure
hold on
for i = 1:length(p)
    pi = p{i};
    ai = a{i};
    if isClose(i) == 0
        C = G2Curve(pi,ai,u,'unclose');
    else
        C = G2Curve(pi,ai);      
    end
    C.drawCurve;
    C.showControlPoints;
end
axis equal; axis off;