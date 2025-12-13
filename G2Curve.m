classdef G2Curve
    properties
        controlPoints           %controlpoints
        pointnum                %controlpoint number
        curvaturePar            %curvature parameter
        curvePar                %0:0.01:1
        kappa                   %directed curvature
        isClose                 
        kappaLength             %curvature scaling parameter(for drawing)
        isFlat                  
        curve                   
        curvature               
        FPoints                 
        normal                  
        tangent                 
        maxdistance             
        blendParameter
        weight                  
        interpPar               
    end
    
    methods
        function obj = G2Curve(p, a, u, isClose, w)
            if nargin < 1
                p = [0.5 0; 1 0.5; 0.5 1; 0 0.5];
            end
            if nargin < 2
                a = ones(length(p), 1);
            end
            if length(a) < numel(p(:,1))
                a = [a;ones(numel(p(:,1))-length(a),1)];
            end
            if length(a) > numel(p(:,1))
                a = a(1:numel(p(:,1)),1);
            end
            if nargin < 3
                u = 0:0.01:1;
            end
            if nargin < 4
                isClose = 'close';
            end
            if nargin < 5
                w = 0.83;
            end

            obj.controlPoints = p;
            if numel(obj.controlPoints(1,:)) == 2
                obj.controlPoints = obj.controlPoints*[1;1i];
            end
            obj.curvaturePar = a;
            obj.curvePar = u;
            obj.isClose = isClose;
            obj.weight = w;

            [obj.curve,obj.kappa,obj.FPoints,obj.maxdistance,obj.normal,obj.tangent,obj.blendParameter,obj.interpPar] = G2(p, u, a, isClose, w);
            obj.curvature = abs(obj.kappa);
            obj.pointnum = length(obj.controlPoints);
            if strcmp(obj.isClose,'unclose')
                obj.pointnum = obj.pointnum-2;
            end
        end 

        function drawCurve(obj,Color,LineWidth)
            if nargin < 3
                LineWidth = 2;
            end
            if nargin < 2
                Color = 'b';
            end
            plot(real(obj.curve),imag(obj.curve),"LineWidth",LineWidth,"Color",Color);
            hold on
        end
        function drawCurvatureCurve(obj,kappalength,track,LineWidth,Color)
            if nargin < 2
                kappalength = 0.3*max(abs(obj.controlPoints-obj.controlPoints([2:length(obj.controlPoints),1])));
            end
            if nargin < 3
                track = 'on';
            end
            if nargin < 4
                LineWidth = 1.2;
            end
            if nargin < 5
                Color = 'm';
            end
            curvaturecurve = obj.curve - kappalength*log10(1+obj.curvature).*obj.kappa./obj.curvature;
            plot(real(curvaturecurve),imag(curvaturecurve),"LineWidth",LineWidth,"Color",Color);
            if strcmp(track,'on')
                d = floor((length(obj.curvePar)-1)/20);
                hold on
                line([real(obj.curve(1:d:end));real(curvaturecurve(1:d:end))],[imag(obj.curve(1:d:end));imag(curvaturecurve(1:d:end))],'Color','g','LineWidth',0.5);
            end
            hold on
        end

        function showControlPoints(obj,Color,isFill,ghostpoints,size)
            if nargin < 2
                Color = 'r';
            end
            if nargin < 3
                isFill = 'unfill';
            end
            if nargin < 4
                ghostpoints = 'off';
            end
            if nargin < 5
                size = [];
            end
            C = obj.controlPoints;
            if strcmp(obj.isClose, 'unclose')
                C = C(2:end-1,:);
                if strcmp(ghostpoints, 'on')
                    scatter(real(obj.controlPoints),imag(obj.controlPoints),[],'g');
                end
            end

            if strcmp(isFill,'filled')
                scatter(real(C),imag(C),size,Color,'filled');
            else
                scatter(real(C),imag(C),size,Color);
            end
        end
        function showFPoints(obj,Color)
            if nargin < 2
                Color = 'm';
            end
            F = obj.FPoints;
            if strcmp(obj.isClose,'unclose')
                F = F(2:end-1,:);
            end
            scatter(real(F(:,1)),imag(F(:,1)),Color,'filled');
            hold on
            scatter(real(F(:,3)),imag(F(:,3)),Color,'filled');
        end

    end
end
function [q,kappa,c,d,n,t,blend,s] = G2(p, u, alpha, isClose, w, h)

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Function for calculating G2-LMC curve
% p: control points, u: parameter, h: curve handle for callback
% varargin contains options: 1. closed curve or not 2. curvature parameter of interpolation function
% p will be converted to complex numbers for calculation

% Initialize curve points q and curvature kappa, used to get F control points c of interpolation function
kappa = [];
c = [];
d = [];
if numel(p(1,:)) == 2
    p = p*[1;1i];
end
q = [];

% Solve for F points and parameters of blending functions
c0 = p([length(p),1:length(p)-1]');
c2 = p([2:length(p),1]');
[c1,s] = getC1S(p,c0,c2);
[n, t] = getNT(c0, c1, c2, s);%[n, t, id] = getNT(c0, c1, c2, s);

id = getFlat(p); % Flat points
id_0 = mod(id-2,length(p))+1;
id_2 = mod(id,length(p))+1;
alpha(intersect(find(alpha<1/0.83),id_0)) = 1/0.83;
alpha(intersect(find(alpha<1/0.83),id_2)) = 1/0.83;
[c0,c2] = getc0c2(p,n,t,alpha,w);
c0(id) = (1-w)*p(id)+w*p(id_0);
c2(id) = (1-w)*p(id)+w*p(id_2);
c2(id_0) = 1./alpha(id_0).*p(id)+(1-1./alpha(id_0)).*p(id_0);
c0(id_2) = 1./alpha(id_2).*p(id)+(1-1./alpha(id_2)).*p(id_2);
[c1,s] = getC1S(p,c0,c2);
[n, t] = getNT(c0, c1, c2, s);
[a, b] = getab(p,c1,c0,c2, n, t, s);
a(id) = 0.5;
b(id) = 0.5;
blend = [a,b];
% According to the option value, choose to draw closed or open curve
if strcmp(isClose, 'close')
    for i = 1:length(p)
        F1 = QuadraticBezier(c0(i),c1(i),c2(i));
        j = mod(i,length(p))+1; % Next point after i
        F2 = QuadraticBezier(c0(j),c1(j),c2(j));
        [xx,yy] = Bezierblend(a(i),b(j));
        C = @(t)F1(s(i)+(1-s(i)).*xx(t)).*yy(t) + F2(s(j).*xx(t)).*(1-yy(t));
        d1 = max(abs(imag(C(u)-p(i)).*real(p(i)-p(j))-imag(p(i)-p(j)).*real(C(u)-p(i))))/abs(p(i)-p(j))^2;
        d = [d;d1];
        q = [q, C(u)];
        kappa = [kappa, curvature(C,u)];   
    end
    c=[c0,c1,c2];    
elseif strcmp(isClose, 'unclose')
    for i = 2:length(p)-2
        F1 = QuadraticBezier(c0(i),c1(i),c2(i));
        F2 = QuadraticBezier(c0(i+1),c1(i+1),c2(i+1));
        [xx,yy] = Bezierblend(a(i),b(i+1));
        C = @(t)F1(s(i)+(1-s(i)).*xx(t)).*yy(t) + F2(s(i+1).*xx(t)).*(1-yy(t));
        q = [q, C(u)];
        kappa = [kappa, curvature(C,u)];
    end
    c=[c0,c1,c2];
else
    error('Unknown style option');
end
if(nargin > 5)
    set(h, 'xdata', real(q), 'ydata',imag(q));
end
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Embedded functions
% Bezierblend(a,b): gets blending function with parameters a,b
% solveQuadpoly(a,b,c): returns the larger positive root of upward-opening quadratic function
% solveQuadpolyN(a,b,c): returns the larger negative root of quadratic function, uses -inf if none exists
% solveQuadpolyP(a,b,c): returns the smaller positive root of quadratic function, uses inf if none exists
% getC1S(p,c0,c2): gets middle F point of interpolation function F and corresponding parameter at control points
% getc0c2(p,n,t,a): gets F points of interpolation function with curvature parameter a, determined by vertex and symmetry axis direction
% getNT(c0, c1, c2, s): gets first and second derivative values at control points of interpolation function
% getab(p,c1,c0,c2, n, t, s): gets parameters of blending function
% QuadraticBezier(p1, p2, p3): gets quadratic Bezier spline
% curvature(F, t): gets directed curvature of curve, represented by complex numbers

function [xx,yy] = Bezierblend(a,b)   
    yy = @(t)(1-t).^4+4*(1-t).^3.*t+3*(1-t).^2.*t.^2;
    xx = @(t)4*a*(1-t).^3.*t+3*(1-t).^2.*t.^2+4.*(1-b).*(1-t).*t.^3+t.^4;   
end
function s = solveQuadpoly(a,b,c)
% f = ax^2+bx+c, a>0, b>0, c<0, s>0
s = (-b+sqrt(b.^2-4.*a.*c))./2./a;
s(imag(s)~=0)=2;
s(isnan(s)) = 2;
end
function s = solveQuadpolyN(a,b,c)
s1 = -2*c./(b+sqrt(b.^2-4.*a.*c));
s2 = -2*c./(b-sqrt(b.^2-4.*a.*c));
s1(imag(s1)~=0) = -inf;
s2(imag(s2)~=0) = -inf;
s1(s1>=0) = -inf;
s2(s2>=0) = -inf;
s=max(s1,s2);
end
function s = solveQuadpolyP(a,b,c)
s1 = -2*c./(b+sqrt(b.^2-4.*a.*c));
s2 = -2*c./(b-sqrt(b.^2-4.*a.*c));
s1(imag(s1)~=0) = inf;
s2(imag(s2)~=0) = inf;
s1(s1<=0) = inf;
s2(s2<=0) = inf;
s=min(s1,s2);
end

function [c1, s] = getC1S(p,c0,c2)
P = length(p);
c1 = zeros(P,1);
s = zeros(P,1);
for i =1:P
    p1 = c0(i);
    p2 = p(i);
    p3 = c2(i);
fun = @(t)(real((p3-p1)*(p3-p1)').*t.^3+3.*real((p3-p1)*(p1-p2)').*t.^2+...
    real((3*p1-2*p2-p3)*(p1-p2)').*t-real((p2-p1)*(p2-p1)'));
x0 = 0.5;
t = fzero(fun, x0);
if t > 1||t < 0
    error('sth wrong')
end
% disp(fun(t));
c1(i) = (p2-(1-t)^2.*p1-t^2.*p3)/(2*(1-t)*t);
s(i) = t;
end

end

function [c0,c2,k0,k2] = getc0c2(p,n,t,a,w)
k = length(p);
alpha = 1-w;
k0c = solveQuadpolyN(real(.5.*a.*n.*conj(t([k,1:k-1]))),...
    real(t.*conj(t([k,1:k-1]))),real((p-p([k,1:k-1])).*conj(t([k,1:k-1]))));
k0d = solveQuadpolyN(real(.5.*a.*n.*conj(n([k,1:k-1]))),...
    real(t.*conj(n([k,1:k-1]))),real((p-p([k,1:k-1])).*conj(n([k,1:k-1]))));
k2c = solveQuadpolyP(real(.5.*a.*n.*conj(t([2:k,1]))),...
    real(t.*conj(t([2:k,1]))),real((p-p([2:k,1])).*conj(t([2:k,1]))));
k2d = solveQuadpolyP(real(.5.*a.*n.*conj(n([2:k,1]))),...
    real(t.*conj(n([2:k,1]))),real((p-p([2:k,1])).*conj(n([2:k,1]))));
k0 = (1-alpha)*max(k0c,k0d);
k2 = (1-alpha)*min(k2c,k2d);
c0 = p+k0.*t+.5*k0.^2.*a.*n;
c2 = p+k2.*t+.5*k2.^2.*a.*n;
for i=1:k
    if n(i)==0
        c0(i) = p(i)*(1-alpha)+p(mod(i-2,k)+1)*alpha;
        c2(i) = p(i)*(1-alpha)+p(mod(i,k)+1)*alpha;
    end
end
end
function [n, t, id] = getNT(c0, c1, c2, s)
n = 2*c0 + 2*c2 - 4*c1;
t = 2*(1-s).*(c1-c0) + 2*s.*(c2-c1);
id=find(abs(n)./abs(c2-c0)<10^-10);
end

function [a,b] = getab(p,c1,c0,c2, n, t, s)
k = length(c1);
e = 0 * 10^(-10);
b0 = real((c1([2:k,1]')-p)./n);    
b2 = real((c1([k,1:k-1]')-p)./n); 
d0 = real((c0([2:k,1]')-p)./n);    
d2 = real((c2([k,1:k-1]')-p)./n); 

a0 = real((c0([2:k,1]')-p).*conj(t))./(t.*conj(t)+e);    
a2 = real((c2([k,1:k-1]')-p).*conj(t))./(t.*conj(t)+e);
a = 1./solveQuadpoly(3/8./s([2:k,1]).*d0.*(1+a0./(1-s)),...
    -4.*d0./3./s([2:k,1]),...
    (1-s)./s([2:k,1]).*a0+2*(d0-b0)-(1-s)./s([2:k,1]));
b = 1./solveQuadpoly(3/8./(1-s([k,1:k-1])).*d2.*(1-a2./s),...
    -4.*d2./3./(1-s([k,1:k-1])),...
    -s./(1-s([k,1:k-1])).*a2+2*(d2-b2)-s./(1-s([k,1:k-1])));
a(a>0.5) = 0.5;
b(b>0.5) = 0.5;
for i = 1:k
    j = mod(i-2,k)+1;
    if d0(i)/(1-s(i))^2/a(i)^2>d2(i)/s(i)^2/b(i)^2
        b(i) = (1-s(i))/s(i)*sqrt(d2(i)/d0(i))*a(i);
    else
        a(i) =s(i)/(1-s(i))*sqrt(d0(i)/d2(i))*b(i);
    end
end
end 

function id = getFlat(p)
p0 = p([length(p),1:length(p)-1])-p;
p2 = p([2:length(p),1])-p;
id = find(abs(real(p0).*real(p2)+imag(p0).*imag(p2))./abs(p0.*p2)>0.99995);
end

function F = QuadraticBezier(p1, p2, p3)
F = @(x) (1-x).^2 * p1 + 2*(1-x) .* x * p2 + x.^2 * p3;
end

function kappa = curvature(F, t)
d = t(2) - t(1);
dft = (F(t+d/2)-F(t-d/2)) / d;
dftt = (F(t+d)-2*F(t)+F(t-d)) / d^2;
kappa = imag(conj(dft).*dftt)./abs(dft).^4.*dft*1i;
end