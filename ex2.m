clear all;
close all

xx = [];
yy = [];
zz = [];
x = [];
y = [];
z = [];

tagdata = load('X:\MT\Mapfilterstuff\output\pointcloud2.txt');

xx = tagdata(:,1);
yy = tagdata(:,2);
zz = tagdata(:,3);
n = 1;
for kk = 1:length(xx)
    if((abs(xx(kk))<5)&&((abs(yy(kk))<100))&&((zz(kk)<5.5)&&(zz(kk)>-1.2)))
%         if ((x(kk)<-0.4)||(abs(y(kk))>23)||(x(kk))>1)
        x(n,1) = xx(kk);
        y(n,1) = yy(kk);
        z(n,1) = zz(kk);
        n = n+1;
%         end
    end

end
% 
% figure(1);
% scatter3(Data(:,1),Data(:,2),Data(:,3));
% 
% Data(n,:) = Data(n,:) + abs(min(Data(n,1)));
% Data(n,:) = Data(n,:) + abs(min(Data(n,2)));

step = 0.05;

xgrid = min(x):step:max(x)+1;
ygrid = min(y):step:max(y)+1;
zgrid = min(z):step:max(z)+1;

xsize  = size(xgrid,2)+1;
ysize  = size(ygrid,2)+1;
zsize  = size(zgrid,2)+1;

xori = min(x);
yori = min(y);
zori = min(z);

Z = zeros(ysize,zsize);

for k=1:length(x)
%     
%     kx = round(Data(k,1)*100) + xsize;
%     ky = round(Data(k,2)*100) + ysize;

    kx = round((y(k,1)-yori)/step)+1;
    ky = round((z(k,1)-zori)/step)+1;
    Z(kx,ky) = Z(kx,ky) +1;
    
end

Z = Z';

% imshow(Z);

theta = 0:0.1:180;
[R,xp] = radon(Z,theta);
imshow(R,[],'Xdata',theta,'Ydata',xp,...
            'InitialMagnification','fit')
xlabel('\theta (degrees)')
ylabel('x''')
colormap(hot), colorbar
iptsetpref('ImshowAxesVisible','off')

[r,c] = find (R == max(max(R)))

rottheta = 90 - c/10
% imshow(imrotate(Z,rottheta))


figure(2);
subplot(2,1,1); imshow(Z);
title('Original Point Cloud')
subplot(2,1,2); imshow(imrotate(Z,rottheta));
title('Rotated Point Cloud')



