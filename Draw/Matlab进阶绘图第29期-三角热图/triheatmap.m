function triheatmap(data,map,xlb,ylb,type)
% data - m*n matrix
% map  - colormap matrix
% xlb     - xticklabels
% ylb     - yticklabels
% type    - 'trid' lower triangle bubble heatmap
%         - 'triu' upper triangle bubble heatmap

% 构造位置数据
x = 1:size(data,2);
y = 1:size(data,1);
[X,Y] = meshgrid(x,y);
Y = flipud(Y);

% 构造网格数据
w = x(2)-x(1);
Xt = (x(1)-w/2):w:(x(end)+w/2);
Yt = (y(1)-y/2):w:(y(end)+y/2);
[Xmesh,Ymesh] = meshgrid(Xt,Yt);
Ymesh = flipud(Ymesh);

% 三角
switch type
    case 'trid'
        % 三角热图数据构造
        temp = triu(ones(size(data)),1);
        idx = find(temp==1);
        data(idx) = NaN;
        Z = data;
        xp = [X(:)-w/2,X(:)+w/2,X(:)+w/2,X(:)-w/2];
        yp = [Y(:)+w/2,Y(:)+w/2,Y(:)-w/2,Y(:)-w/2];
        % 三角网格数据构造
        Zmesh = zeros(size(Xmesh));
        temp = triu(ones(size(Zmesh)),2);
        idx2 = find(temp==1);
        Zmesh(idx2) = NaN;   
        % 颜色构造
        minc = min(data(:));
        maxc = max(data(:));
        num = size(map,1);
        x = linspace(minc,maxc,num);
        cmap = interp1(x,map,data(:),'linear'); %...插值
        % 三角热图绘制
        for i = 1:length(data(:))
            if isnan(data(i))
                patch(xp(i,1:4),yp(i,1:4),[1 1 1],'edgecolor','none')
            else
                patch(xp(i,1:4),yp(i,1:4),cmap(i,1:3),'edgecolor','none')
            end
        end
        hold on
        mesh(Xmesh,Ymesh,Zmesh,'facecolor','none','EdgeColor',[0.2 0.2 0.2],'LineWidth',1)
        view(0,90)
        axis equal tight off
        set(gca,'xticklabels',[])
        for i = 1:size(Z,1)
            text(0.3,Y(i,1),xlb{i},'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',11,'FontName','arial','color','k')
            text(X(1,i),Y(i,1)+0.6,ylb{i},'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',11,'FontName','arial','color','k')
        end
        colorbar('southoutside')
        caxis([minc maxc])

    case 'triu'
        % 三角热图数据构造
        temp = triu(ones(size(data)));
        idx = find(temp==0);
        data(idx) = NaN;
        Z = data;
        xp = [X(:)-w/2,X(:)+w/2,X(:)+w/2,X(:)-w/2];
        yp = [Y(:)+w/2,Y(:)+w/2,Y(:)-w/2,Y(:)-w/2];
        % 三角网格数据构造
        Zmesh = zeros(size(Xmesh));
        temp = triu(ones(size(Zmesh)),-1);
        idx2 = find(temp==0);
        Zmesh(idx2) = NaN;   
        % 颜色构造
        minc = min(data(:));
        maxc = max(data(:));
        num = size(map,1);
        x = linspace(minc,maxc,num);
        cmap = interp1(x,map,data(:),'linear'); %...插值
        % 三角热图绘制
        for i = 1:length(data(:))
            if isnan(data(i))
                patch(xp(i,1:4),yp(i,1:4),[1 1 1],'edgecolor','none')
            else
                patch(xp(i,1:4),yp(i,1:4),cmap(i,1:3),'edgecolor','none')
            end
        end
        hold on
        mesh(Xmesh,Ymesh,Zmesh,'facecolor','none','EdgeColor',[0.2 0.2 0.2],'LineWidth',1)
        view(0,90)
        axis equal tight off
        set(gca,'yticklabels',[])
        for i = 1:size(Z,1)
            text(X(1,i),Y(1)+0.5,xlb{i},'HorizontalAlignment','center','VerticalAlignment','bottom','FontSize',11,'FontName','arial','color','k')
            text(X(1,i)-0.6,Y(i,1),ylb{i},'HorizontalAlignment','right','VerticalAlignment','middle','FontSize',11,'FontName','arial','color','k')
        end
        colorbar('eastoutside')
        caxis([minc maxc])
end

end