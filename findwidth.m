function  [Distance mid xleft xright]= findwidth(cloudpoint)

% Function   : findwidth
% 
% Purpose    : Calculate width based on Randon transform
% 
% Parameters : cloudpoint, Point cloud data
% 
% Return     : Distance, Corridor width
%              mid, Center of the corridor
%              xleft, xrightL left and right corridor



step = 0.05;

x(:,1) = cloudpoint(:,1);
y(:,1) = cloudpoint(:,2);
z(:,1) = cloudpoint(:,3);

xgrid = min(x):step:max(x)+1;
ygrid = min(y):step:max(y)+1;
zgrid = min(z):step:max(z)+1;

xsize  = size(xgrid,2)+1;
ysize  = size(ygrid,2)+1;
zsize  = size(zgrid,2)+1;

xori = min(x);
yori = min(y);
zori = min(z);

Z = zeros(xsize,ysize);



for k=1:length(x)
%     
%     kx = round(Data(k,1)*100) + xsize;
%     ky = round(Data(k,2)*100) + ysize;

    kx = round((x(k,1)-xori)/step)+1;
    ky = round((y(k,1)-yori)/step)+1;
    Z(kx,ky) = Z(kx,ky) +1;
    
end

theta = 0:0.1:180;
[R,xp] = radon(Z,theta);

figure;title('Radon Transform')
imshow(R,[],'Xdata',theta,'Ydata',xp,...
            'InitialMagnification','fit')
xlabel('\theta (degrees)')
ylabel('x''')
colormap(hot), colorbar
iptsetpref('ImshowAxesVisible','off')

[r,c] = find (R == max(max(R)));

c1 = c;
% c2 = c1 + 900
% figure(33);
% plot(R(:,c1))
% figure(44);
% plot(R(:,c2))

varargout = peakfinder(R(:,c1),max(R(:,c1))/6, 1)

Max = max(varargout);

sec = max(varargout(varargout~=max(varargout)));

Distance = abs((Max - sec)*step)

mid = abs((Max - sec)/2) + sec;

midpos =  min(cloudpoint(:,1)) + (mid/length(R(:,c1)))*(max(cloudpoint(:,1))-min(cloudpoint(:,1)));

xleft = midpos - Distance/2
xright = midpos + Distance/2
end