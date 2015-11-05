function [posBuff,slope,xIntercept,prevLines] = chooseXY(lines,posBuff,frameCenter,x0,y0,slope,xIntercept,Hough,lane,prevLines,wRho,wTheta,AverageWeight)

n = struct([]);

if strcmp(lane,'left')
    j = 1;
    for i = 1:numel(lines)
        p = lines(i).point1;
        if p(1) < frameCenter(2)
            n(j).point1 = lines(i).point1;
            n(j).point2 = lines(i).point2;
            n(j).rho = lines(i).rho;
            n(j).theta = lines(i).theta;
            j = j+1;
        end
    end
    lines = n;
else
    j = 1;
    for i = 1:numel(lines)
        p = lines(i).point2;
        if  p(1) > frameCenter(2)
            n(j).point1 = lines(i).point1;
            n(j).point2 = lines(i).point2;
            n(j).rho = lines(i).rho;
            n(j).theta = lines(i).theta;
            j = j+1;
        end
    end
    lines = n;
end;

if numel(lines) == 0
    if posBuff(2,end) < 3
        lines = prevLines;
        posBuff(2,end) = posBuff(2,end) + 1;
    else
        prevLines = struct([]);
    end;
else
    posBuff(2,end) = 0;
    
    xypoints1 = reshape([lines.point1],[2,length(lines)]);
    xypoints2 = reshape([lines.point2],[2,length(lines)]);
    
    rhos = reshape([lines.rho],[1,length(lines)]);
    thetas = reshape([lines.theta],[1,length(lines)]);
    
    prevD = posBuff(1,end);
    prevRho = posBuff(3,end);
    prevTheta = posBuff(4,end);
    
    vRho = abs(rhos - prevRho);
    vTheta = abs(thetas - prevTheta);
    
    posBuff(:,1:end-1) = posBuff(:,2:end);
    
    HY = max(uint16(size(Hough,1)/2+prevRho-wRho),1):min(uint16(size(Hough,1)/2+prevRho+wRho),size(Hough,1));
    HX = max(uint16(size(Hough,2)/2+prevTheta-wTheta),1):min(uint16(size(Hough,2)/2+prevTheta+wTheta),size(Hough,2));
        
    if sum(vRho < 20 & vTheta < 4) > numel(vRho)

    elseif sum(sum(Hough(HY,HX)))/((2*(wRho-15/prevD)+1)*(2*(wTheta-15/prevD)+1)) > AverageWeight
        
    else
        
        X1 = xypoints1(1,:);
        Y1 = xypoints1(2,:);
        X2 = xypoints2(1,:);
        Y2 = xypoints2(2,:);
        
        Q1 = [X1',Y1',zeros(size(X1,2),1)];
        Q2 = [X2',Y2',zeros(size(X2,2),1)];
        P = repmat([x0 y0 0],size(X2,2),1);
        a = Q1-Q2;
        b = P-Q2;
        for i = 1:size(P,1)
            d(i) = norm(cross(a(i,:),b(i,:)))/norm(a(i,:));
        end
        
        [posBuff(1,end),k] = min(d);
        
        
        xy = [lines(k).point1; lines(k).point2];
        
        slope = (xy(2,1)-xy(1,1))/(xy(2,2)-xy(1,2));
        
        xIntercept = xy(1,1) - slope*xy(1,2);
        
        posBuff(3,end) = lines(k).rho;
        posBuff(4,end) = lines(k).theta;
    end;
end;

prevLines = lines;

end

% function [posBuff,slope,xIntercept,prevLines] = chooseXY(lines,posBuff,frameCenter,x0,y0,slope,xIntercept,Hough,lane,prevLines)
% 
% wRho = 3;
% wTheta = 3;
% 
% n = struct([]);
% 
% if strcmp(lane,'left')
%     j = 1;
%     for i = 1:numel(lines)
%         p = lines(i).point1;
%         if p(1) < frameCenter(2)
%             n(j).point1 = lines(i).point1;
%             n(j).point2 = lines(i).point2;
%             n(j).rho = lines(i).rho;
%             n(j).theta = lines(i).theta;
%             j = j+1;
%         end
%     end
%     lines = n;
% else
%     j = 1;
%     for i = 1:numel(lines)
%         p = lines(i).point2;
%         if  p(1) > frameCenter(2)
%             n(j).point1 = lines(i).point1;
%             n(j).point2 = lines(i).point2;
%             n(j).rho = lines(i).rho;
%             n(j).theta = lines(i).theta;
%             j = j+1;
%         end
%     end
%     lines = n;
% end;
% 
% if numel(lines) == 0
%     if posBuff(2,end) < 3
%         lines = prevLines;
%         posBuff(2,end) = posBuff(2,end) + 1;
%     else
%         prevLines = struct([]);
%     end;
% else
%     posBuff(2,end) = 0;
%     
%     xypoints1 = reshape([lines.point1],[2,length(lines)]);
%     xypoints2 = reshape([lines.point2],[2,length(lines)]);
%     
%     rhos = reshape([lines.rho],[1,length(lines)]);
%     thetas = reshape([lines.theta],[1,length(lines)]);
%     
%     prevRho = posBuff(3,end);
%     prevTheta = posBuff(4,end);
%     
%     vRho = abs(rhos - prevRho);
%     vTheta = abs(thetas - prevTheta);
%     
%     posBuff(:,1:end-1) = posBuff(:,2:end);
%     
%     HY = max(uint16(size(Hough,1)/2+prevRho-wRho),1):min(uint16(size(Hough,1)/2+prevRho+wRho),size(Hough,1));
%     HX = max(uint16(size(Hough,2)/2+prevTheta-wTheta),1):min(uint16(size(Hough,2)/2+prevTheta+wTheta),size(Hough,2));
%         
%     if sum(vRho < 20 & vTheta < 4) > numel(vRho)
% 
%     elseif sum(sum(Hough(HY,HX)))/((2*wRho+1)*(2*wTheta+1)) > 10
%         % 
%         
%     else
%         
%         X1 = xypoints1(1,:);
%         Y1 = xypoints1(2,:);
%         X2 = xypoints2(1,:);
%         Y2 = xypoints2(2,:);
%         
%         Q1 = [X1',Y1',zeros(size(X1,2),1)];
%         Q2 = [X2',Y2',zeros(size(X2,2),1)];
%         P = repmat([x0 y0 0],size(X2,2),1);
%         a = Q1-Q2;
%         b = P-Q2;
%         for i = 1:size(P,1)
%             d(i) = norm(cross(a(i,:),b(i,:)))/norm(a(i,:));
%         end
%         
%         [~,k] = min(d);
%         
%         
%         xy = [lines(k).point1; lines(k).point2];
%         
%         slope = (xy(2,1)-xy(1,1))/(xy(2,2)-xy(1,2));
%         
%         xIntercept = xy(1,1) - slope*xy(1,2);
%         
%         posBuff(3,end) = lines(k).rho;
%         posBuff(4,end) = lines(k).theta;
%     end;
% end;
% 
% prevLines = lines;
% 
% end



%
% xypoints1 = reshape([lines.point1],[2,length(lines)]);
% xypoints2 = reshape([lines.point2],[2,length(lines)]);
%
% rhos = reshape([lines.rho],[1,length(lines)]);
% thetas = reshape([lines.theta],[1,length(lines)]);
%
% prevRho = posBuff(3,end);
% prevTheta = posBuff(4,end);
%
% vRho = abs(rhos - prevRho);
% vTheta = abs(thetas - prevTheta);
%
% posBuff(:,1:end-1) = posBuff(:,2:end);
%
% if sum(vRho < 5 & vTheta < 4) > 0
%
% else
%
%     X1 = xypoints1(1,:);
%     Y1 = xypoints1(2,:);
%     X2 = xypoints2(1,:);
%     Y2 = xypoints2(2,:);
%
%     a = (Y2-Y1);
%     b = (X2-X1);
%     c = (Y1.*X2-X1.*Y2);
%     d = abs(a.*x0+b.*y0+c)./sqrt(a.^2+b.^2);
%
%     [~,k] = min(d);
%
%
%     xy = [lines(k).point1; lines(k).point2];
%
%     slope = (xy(2,1)-xy(1,1))/(xy(2,2)-xy(1,2));
%
%     xIntercept = xy(1,1) - slope*xy(1,2);
%
%     posBuff(3,end) = lines(k).rho;
%     posBuff(4,end) = lines(k).theta;
%
% end;
%
% end
%










% function [posBuff,slope,xIntercept] = chooseXY(lines,posBuff,frameCenter,x0,y0,slope,xIntercept)
%
% xypoints=reshape([lines.point1],[2,length(lines)]);
% xypoints = xypoints';
% ypoints = xypoints(:,2);
%
% [~,k] = max(ypoints);
%
% xy = [lines(k).point1; lines(k).point2];
%
% xypoints1 = reshape([lines.point1],[2,length(lines)]);
% X1 = xypoints1(1,:);
% Y1 = xypoints1(2,:);
% xypoints2 = reshape([lines.point2],[2,length(lines)]);
% X2 = xypoints2(1,:);
% Y2 = xypoints2(2,:);
%
% x0 = frameCenter(2);
% y0 = frameCenter(1)*2 - frameCenter(1)/4;
%
% a = (Y2-Y1);
% b = (X2-X1);
% c = (Y1.*X2-X1.*Y2);
% d = abs(a.*x0+b.*y0+c)./sqrt(a.^2+b.^2);
%
% [~,k] = min(d);
%
% if abs(lines(k).rho - posBuff(2,end)) > 5
%
% xy = [lines(k).point1; lines(k).point2];
%
% slope = (xy(2,1)-xy(1,1))/(xy(2,2)-xy(1,2));
%
% xIntercept = xy(1,1) - slope*xy(1,2);
%
% yIntercept = int16(-xIntercept/slope);
%
%
% if ((xy(1,2) - median(posBuff(1,:))) > 10)
%     xy(1,2) = median(posBuff(1,:));
%     xy(1,1) = posBuff(2, find(posBuff(1,:) == xy(1,2), 1));
% end;
%
% posBuff(:,1:end-1) = posBuff(:,2:end);
% posBuff(1,end) = yIntercept;
% posBuff(2,end) = lines(k).rho;
%
% xIntercept = -slope*median(posBuff(1,:));
%
% else
%
%
%
% end
%
% posBuff(:,1:end-1) = posBuff(:,2:end);
% posBuff(1,end) = xy(1,2);
% posBuff(2,end) = xy(1,1);
%
% end