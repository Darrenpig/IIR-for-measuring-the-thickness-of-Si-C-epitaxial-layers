function [s,idx] = BubbleScatter(data,sz,fc)
% data - m*1向量，x坐标
% sz   - m*1向量，尺寸参数
% fc   - 特征渲染滑珠气泡图：m*1向量，颜色参数
%        其它：1*3颜色数据
% Author:AKun
% 公众号：阿昆的科研日常

N = length(data);
y = 1:N;

[Sdata,idx] = sort(data);

if length(sz) == 1
    s = scatter(Sdata, y, sz, 'k', 'filled');
elseif size(fc,1)==1
    sz = sz(:);
    sz = sz(idx);
    s = bubblechart(Sdata, y, sz, fc);
elseif size(fc,1)>1
    sz = sz(:);
    sz = sz(idx);
    fc = fc(:);
    fc = fc(idx); 
    s = bubblechart(Sdata, y, sz, fc);
end

end