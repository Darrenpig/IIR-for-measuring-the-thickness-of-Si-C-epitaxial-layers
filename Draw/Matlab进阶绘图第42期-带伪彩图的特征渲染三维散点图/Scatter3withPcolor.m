function Scatter3withPcolor(data,sz,cellsize,zpos)
% data          -  m*3数据矩阵
% sz            -  散点尺寸
% cellsize      -  网格尺寸
% zpos          -  伪彩图位置
% By 阿昆的科研日常

% 初始化
x = data(:,1);
y = data(:,2);
v = data(:,3);
% 生成伪彩图格网数据
[xmesh,ymesh,vmesh] = kDSM(data,cellsize);
% 绘制
scatter3(x, y, v, sz, v, 'filled')
hold on
surf(xmesh,ymesh,zpos*ones(size(xmesh)),vmesh,'EdgeColor','none');
end

function [xmesh,ymesh,vmesh] = kDSM(data,cellsize)
% data          -  m*3数据矩阵
% cellsize      -  网格尺寸
% By 阿昆的科研日常
x = data(:,1);
y = data(:,2);
v = data(:,3);
minx = min(x);
maxx = max(x);
miny = min(y);
maxy = max(y);
xi = minx:cellsize:maxx;
yi = miny:cellsize:maxy;
[xmesh,ymesh] = meshgrid(xi,yi);
vmesh = griddata(x,y,v,xmesh,ymesh);
end