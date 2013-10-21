


d = load('X:\MT\Mapfilterstuff\output\pointcloud2.txt');


step = 0.05;
xx = [];
yy = [];
zz = [];
x = [];
y = [];
z = [];
X_corr = [];
xx = d(:,1);
yy = d(:,2);
zz = d(:,3);
n=1;
for kk = 1:length(xx)
    if((abs(xx(kk))<5)&&((abs(yy(kk))<100))&&((zz(kk)<5.5)&&(zz(kk)>-1.2)))
%         if ((x(kk)<-0.4)||(abs(y(kk))>23)||(x(kk))>1)
        x(n,1) = xx(kk);
        y(n,1) = yy(kk);
        z(n,1) = zz(kk);
        data(n,:) = [x(n,1) y(n,1) z(n,1)];
        n = n+1;
%         end
    end

end


figure(55);
scatter3(data(:,1),data(:,2),data(:,3));
title('3D Point Cloud');
xlabel('X (m)');
ylabel('Y (m)');


%  step = 0.05;

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

c1 = c
c2 = c1 + 900
figure(3333);
plot(R(:,c1))
figure(4444);
plot(R(:,c2))

rotthetaxy = 90-c/10
thetaxy = rotthetaxy;

% figure;
% imshow(Z);
% % title('Original Point Cloud')
% % subplot(2,1,2); imshow(imrotate(Z,thetaxy));
% % title('Rotated Point Cloud')
figure;
imshow(imrotate(Z,thetaxy));

Rxy = [cosd(thetaxy) -sind(thetaxy); sind(thetaxy) cosd(thetaxy)];
 
  for k=1:length(data(:,1))
    
    vv = [data(k,2) data(k,1)] * Rxy;
    data(k,2) = vv(1,1);
    data(k,1) = vv(1,2);
%         plot3(x(k), y(k), z(k), 'ws--', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
  end
 
  nnn1 = 1;
  nnn2 = 1;
   for k=1:length(data(:,1))
    
        if( data(k,1) <= 0)
               datal(nnn1,2) = data(k,2);
               datal(nnn1,1) = data(k,1);
               nnn1 = nnn1 +1;
        else
               datar(nnn2,2) = data(k,2);
               datar(nnn2,1) = data(k,1);
               nnn2 = nnn2 +1;
        end
%         plot3(x(k), y(k), z(k), 'ws--', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
   end
   

Z1 = zeros(xsize,ysize);

for k=1:length(datar)
%     
%     kx = round(Data(k,1)*100) + xsize;
%     ky = round(Data(k,2)*100) + ysize;

    kx = round((datar(k,1)-xori)/step)+1;
    ky = round((datar(k,2)-yori)/step)+1
    Z1(kx,ky) = Z1(kx,ky) +1;
    
end
figure;
imshow(Z);


theta = 0:0.1:180;
[R,xp] = radon(Z1,theta); 

figure;title('Radon Transform')
imshow(R,[],'Xdata',theta,'Ydata',xp,...
            'InitialMagnification','fit')
xlabel('\theta (degrees)')
ylabel('x''')
colormap(hot), colorbar
iptsetpref('ImshowAxesVisible','off')

[r,c] = find (R == max(max(R)));

c1 = c
c2 = c1 + 900
figure(33);
plot(R(:,c2))
figure(66);
plot(R(:,c1))
figure(663);
plot(datar(:,1),datar(:,2))



Z1 = zeros(xsize,ysize);
close all;
for k=1:length(datal)
%     
%     kx = round(Data(k,1)*100) + xsize;
%     ky = round(Data(k,2)*100) + ysize;

    kx = round((datal(k,1)-xori)/step)+1;
    ky = round((datal(k,2)-yori)/step)+1
    Z1(kx,ky) = Z1(kx,ky) +1;
    
end

theta = 0:0.1:180;
[R,xp] = radon(Z1,theta); 

figure;title('Radon Transform')
imshow(R,[],'Xdata',theta,'Ydata',xp,...
            'InitialMagnification','fit')
xlabel('\theta (degrees)')
ylabel('x''')
colormap(hot), colorbar
iptsetpref('ImshowAxesVisible','off')

[r,c] = find (R == max(max(R)));

c1 = c
c2 = c1 + 900
figure(33);
plot(R(:,c2))
figure(66);
plot(R(:,c1))
figure(663);
plot(datar(:,1),datar(:,2))
