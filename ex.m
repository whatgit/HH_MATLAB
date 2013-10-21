clear all;
close all

xx = [];
yy = [];
zz = [];
x = [];
y = [];
z = [];

CorrWidth = 2.90;
PillarWidth = 2.70;

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

Z = zeros(xsize,ysize);

for k=1:length(x)
%     
%     kx = round(Data(k,1)*100) + xsize;
%     ky = round(Data(k,2)*100) + ysize;

    kx = round((x(k,1)-xori)/step)+1;
    ky = round((y(k,1)-yori)/step)+1;
    Z(kx,ky) = Z(kx,ky) +1;
    
end

% Z = Z';
 
imshow(Z);
%%
theta = 0:0.1:180;
[R,xp] = radon(Z,theta);

figure(1);title('Radon Transform')
imshow(R,[],'Xdata',theta,'Ydata',xp,...
            'InitialMagnification','fit')
xlabel('\theta (degrees)')
ylabel('x''')
colormap(hot), colorbar
iptsetpref('ImshowAxesVisible','off')

[r,c] = find (R == max(max(R)))

c1 = c
c2 = c1 + 900
% figure(33);
% subplot(2,2,[1 3]); imshow(Z);
% title('Point Cloud Porjection on x-y plane');
% xlabel('X');
% ylabel('Y');
% 
% 
% subplot(2,2,2);plot(R(:,c1))
% title('Randon Transform');
% xlabel('Pixel');
% ylabel('Intensity');
% 
% 
% subplot(2,2,4);plot(R(:,c2))
% title('Randon Transform 90 degree');
% xlabel('Pixel');
% ylabel('Intensity');



% 
% nn = 1;
% 
% thresh = max(R(:,c1))/2.5;
% for kk = 1:length(R(:,c1))
%      if(R(kk,c1)>=thresh)
%          NewR(kk,1) = R(kk,c1);
%      else
%          NewR(kk,1) = 0;
%      end   
% end
% 
% varargout = peakfinder(x0, sel, thresh, extrema)
% 
% figure(556);
% plot(NewR(kk,1))

rottheta = -c/10




% imshow(imrotate(Z,rottheta))

%%
figure(2);
subplot(1,2,1); imshow(Z);
title('Original Point Cloud')
subplot(1,2,2); imshow(imrotate(Z,rottheta));
title('Rotated Point Cloud')

%%

RotZ = imrotate(Z,rottheta);

CorrWidth = 2.90;
PillarWidth = 2.70;
step = 0.05;
% Pnum = length(rightp,1);
Pnum = 12;

xgrid = -0.2-CorrWidth/2:step:CorrWidth+0.2;
ygrid = -0.2-PillarWidth*Pnum:step:PillarWidth*Pnum+0.2;



mu = [round(CorrWidth/(step*2)) round(CorrWidth/(step*2))];
% Sigma = [.25 .3; .3 1];
xmesh = 0:step:CorrWidth+0.4;
ymesh = 0:step:CorrWidth+0.4;

[X1,X2] = meshgrid(xmesh,ymesh);
F = mvnpdf([X1(:) X2(:)],mu);
F = reshape(F,length(ymesh),length(xmesh));
figure(3);
surf(xmesh,ymesh,F);
caxis([min(F(:))-.5*range(F(:)),max(F(:))]);
axis([-3 3 -3 3 0 .4])
xlabel('x1'); ylabel('x2'); zlabel('Probability Density');



