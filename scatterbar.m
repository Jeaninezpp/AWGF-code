function scatterbar(x,y,z,scale)
%   根据散点数据绘制3维彩色柱状图
%   scatterbar(x,y,z,scale)  x,y,z是实值数组，用来指定柱子顶面中心点三维坐标。
%              scale是大于0的标量，用来指定柱子的粗细，scale越大，柱子越细，默
%              认情况下根据坐标点自动计算柱子的粗细。
%
%   CopyRight:xiezhh（谢中华）
%   2011.10.31
%   Example：
%       [x,y] = meshgrid(-6:6,-3:0.5:3);
%       z = mvnpdf([x(:),y(:)],[0,0],[4,0;0,1]);
%       scatterbar(x,y,z)
%       scatterbar(x,y,z,50);

% 输入参数类型判断
if nargin < 3
    error('至少需要三个输入参数');
end
if ~isreal(x) || ~isreal(y) || ~isreal(z)
    error('前三个输入应为实值数组');
end

% 提取x,y,z等长部分的元素
x = x(:);
y = y(:);
z = z(:);
n = min([numel(x) numel(y) numel(z)]);
x = x(1:n);
y = y(1:n);
z = z(1:n);

% 计算极差和差分值
rx = range(x);
ry = range(y);
dx = abs(diff(x));
dx = min(dx(dx>0));
dy = abs(diff(y));
dy = min(dy(dy>0));

% 自动计算柱子的粗细
if nargin == 3
    if ~isempty(dx)
        hx = dx/2;
    else
        hx = 0.5;
    end
    if ~isempty(dy)
        hy = dy/2;
    else
        hy = 0.5;
    end
end

% 根据用户输入参数scale计算柱子的粗细
if nargin == 4
    if ~isreal(scale) || scale < 0
        error('第四个输入应为正的标量');
    end
    if rx == 0 && ry == 0
        rx = 0.5*scale;
        ry = rx;
    elseif rx == 0 || ry == 0
        rx = max(rx,ry);
        ry = rx;
    end
    hx = rx/scale;
    hy = ry/scale;
end

% 通过循环绘制三维彩色柱状图
figure
hold on
Xp = [];
Yp = [];
Zp = [];
for i = 1:n
    [xp,yp,zp] = Vertices(x(i),y(i),z(i));
    Xp = [Xp;xp];
    Yp = [Yp;yp];
    Zp = [Zp;zp];
end
%通过surf函数生成彩色的立方体盒子
h = surf(Xp,Yp,Zp,Zp,'FaceColor','interp');
%set(h,'FaceAlpha',0.25);    %设置立方体盒子透明度
grid on
view(3)
hold off

%--------------------------------------------------
% 求柱子顶点的子函数
%--------------------------------------------------
function [xp,yp,zp] = Vertices(x,y,z)
    % 由长方体底面中心坐标求顶点坐标
    xp = [x-hx x-hx x+hx x+hx x-hx
          x-hx x-hx x+hx x+hx x-hx
          x-hx x-hx x+hx x+hx x-hx
          x-hx x-hx x+hx x+hx x-hx
          x    x    x    x    x
          NaN  NaN  NaN  NaN  NaN];
    yp = [y-hy y+hy y+hy y-hy y-hy
          y-hy y+hy y+hy y-hy y-hy
          y-hy y+hy y+hy y-hy y-hy
          y-hy y+hy y+hy y-hy y-hy
          y    y    y    y    y
          NaN  NaN  NaN  NaN  NaN];
    zp = [repmat(linspace(0,z,4)',[1,5]);z z z z z;NaN NaN NaN NaN NaN];
end
end