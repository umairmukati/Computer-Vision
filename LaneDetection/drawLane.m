function image = drawLane(slope, xIntercept, image,lane,mSize)

IMSize = uint16(size(image));

if strcmp(lane,'left')
    yIntercept = -xIntercept/slope;
    if (yIntercept*2) <= IMSize(1)
        x = 1;
        y = 2*yIntercept;
    else
        y = double(IMSize(1));
        x = 2*xIntercept + y*slope;
    end
    slope = 1/slope;
    
    dir = atand(slope);
    
    x = floor(x + [0:0.5:mSize].*cosd(dir));
    y = floor(y + [0:0.5:mSize].*sind(dir));
else
    y = (double(IMSize(2)) - xIntercept*2) / slope;
    if y <= IMSize(1)
        x = double(IMSize(2));
    else
        y = double(IMSize(1));
        x = 2*xIntercept + y*slope;
    end;
    slope = 1/slope;
    
    dir = atand(slope);
    
    x = floor(x - [0:0.5:mSize].*cosd(dir));
    y = floor(y - [0:0.5:mSize].*sind(dir));
end;

% B = s/sqrt(1+slope^2);
% y2 = [y+B,y-B];
% x2 = x+slope*(y2-y)


H = [0 255 0];

for i = 1:size(x,2)
    image(max(y(i)-3,1):min(y(i)+3,IMSize(1)),max(x(i)-3,1):min(x(i)+3,IMSize(2)),:) = ...
        permute(repmat(H,[size(max(y(i)-3,1):min(y(i)+3,IMSize(1)),2),1,size(max(x(i)-3,1):min(x(i)+3,IMSize(2)),2)]),[1,3,2]);
end;

end

% IMSize = uint16(size(image));
%
% Yline = 1:0.5:double(IMSize(1));
% Xline = slope*Yline + xIntercept*frameProcFactor;
%
% [~,b] = find(Xline >= 1 & Xline <= IMSize(2));
% Yline = uint16(round(Yline(b)));
% Xline = uint16(round(Xline(b)));
%
% IM1 = zeros(IMSize(1,1:2));
% for i = 1:size(Yline,2)
%     IM1(max(Yline(i)-3,1):min(Yline(i)+3,IMSize(1)),max(Xline(i)-3,1):min(Xline(i)+3,IMSize(2))) = 1;
% end;
%
% IM1 = IM1.*mask;
%
% image(:,:,1) = image(:,:,1).*uint8((1-IM1));
% image(:,:,2) = image(:,:,2).*uint8((1-IM1));
% image(:,:,3) = image(:,:,3).*uint8((1-IM1));
% image(:,:,2) = image(:,:,2) + uint8(255*IM1);
%
% end