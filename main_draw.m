clc
clear
load('.\data\data_elephant.mat');
%% draw 

u = 0:0.01:1;
figure
for i = 1:length(p)
if isClose(i) == 1
    isclose_cur = 'close';
else
    isclose_cur = 'unclose';
end
C = G2Curve(p{i},a{i},u,isclose_cur);
C.drawCurve;
C.showColorControlPoints;
end
axis equal; axis off

