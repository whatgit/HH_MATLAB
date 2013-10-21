PillarWidth = 2.95;
PillarLength = 2.70;

 scaleWidth = 0.9492;
 scaleLength = 0.9489;
 
 FixedW = PillarWidth*scaleWidth;
 FixedL = PillarLength*scaleLength;


%% Mesh Gird Generate
clear all;

stepm = 0.02;
mu = [0 0];
Sigma = [.005 0;0 0.01]
% Sigma = [1 .5; .5 1];
x1 = -0.2:stepm:.2; x2 = -0.2:stepm:.2;
[X1,X2] = meshgrid(x1,x2);
F = mvnpdf([X1(:) X2(:)],mu,Sigma);
F = reshape(F,length(x2),length(x1));
surf(x1,x2,F);
caxis([min(F(:))-.5*range(F(:)),max(F(:))]);
axis equal
xlabel('x1'); ylabel('x2'); zlabel('Probability Density');

bigmeshx = -1.65:stepm:1.65;
bigmeshy = -0.21:stepm:32.21;
[X1,X2] = meshgrid(x1,x2);

FF = zeros(length(bigmeshx)+1,length(bigmeshy)+1);
  

leftN = 11;
RightN = 11;
lengthF = length(F);

% every Grid 0.01m
CenterSy = 0;
for kk = 1:leftN
    CenterSx = 0;
    
    for k1 = 1:lengthF
        for k2 = 1:lengthF
        FF(CenterSx+k1,CenterSy+k2) = F(k1,k2);
        end
    end
    CenterSy = round(CenterSy + 2.75/stepm);
    kk
end

 CenterSy = 0;
for kk = 1:leftN
    CenterSx = round(2.9/stepm);
       for k1 = 1:lengthF
        for k2 = 1:lengthF
        FF(CenterSx+k1,CenterSy+k2) = F(k1,k2);
        end
       end
    CenterSy = round(CenterSy + 2.75/stepm);
    kk
end

figure;
imshow(FF);

i = 1;
% FF = reshape(FF,length(bigmeshx),length(bigmeshy));
% surf(bigmeshx,bigmeshy,FF);
% caxis([min(F(:))-.5*range(F(:)),max(F(:))]);
% axis equal
% xlabel('x1'); ylabel('x2'); zlabel('Probability Density');



%% Gridlize point cloud

d = load('X:\MT\Mapfilterstuff\output\pointcloud2.txt'); 

step = 0.05;
[data r c] = transOrth(d,step); % Rotate Data till Dominant Direction Arg: points,step

        x(:,1) = data(:,1);
        y(:,1) = data(:,2);
        z(:,1) = data(:,3);

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

    kx = round((data(k,1)-xori)/step)+1;
    ky = round((data(k,2)-yori)/step)+1;
    Z(kx,ky) = Z(kx,ky) +1;
    
end

figure;
imshow(Z);

cur = zeros(size(Z,1)-size(FF,1),(size(Z,2)-size(FF,2)));

for shiftrow = 1: (size(Z,1)-size(FF,1))
    gridxidx = shiftrow
    for shiftcol = 1: (size(Z,2)-size(FF,2))
       gridyidy = shiftcol;
       curin = zeros(size(FF,1),size(FF,2));
       for shiftinrow = 1:size(FF,1)
           for shiftincol = 1:size(FF,2)
                   curin(shiftinrow,shiftincol) = Z(shiftrow+shiftinrow,shiftcol+shiftincol).*FF(shiftinrow,shiftincol);
                   %    Mul
           end
       end
       cur(shiftrow,shiftcol) = sum(sum(curin));
    end
end

[mr,mc] = find (cur == max(max(cur)));

Zres = zeros(size(Z,1),size(Z,2));
    resxidx = mr;
    resyidx = mc;
for kk = 1:size(FF,1)
    for ik = 1:size(FF,2)
        Zres(resxidx+kk,resyidx+ik) = FF(kk,ik);
        
    end
end 

figure;
subplot(3,1,1);imshow(Z);
subplot(3,1,2);imshow(FF);
subplot(3,1,3);imshow(Zres);

%% cal Result

leftpx = min(x) + resxidx*stepm;
leftpy = min(y) + resyidx*stepm;



% 
% step = 0.01;
% xx = [];
% yy = [];
% zz = [];
% x = [];
% y = [];
% z = [];
% X_corr = [];
% xx = d(:,1);
% yy = d(:,2);
% zz = d(:,3);
% n=1;
% for kk = 1:length(xx)
%     if((abs(xx(kk))<5)&&((abs(yy(kk))<100))&&((zz(kk)<5.5)&&(zz(kk)>-1.2)))
% %         if ((x(kk)<-0.4)||(abs(y(kk))>23)||(x(kk))>1)
%         x(n,1) = xx(kk);
%         y(n,1) = yy(kk);
%         z(n,1) = zz(kk);
%         data(n,:) = [x(n,1) y(n,1) z(n,1)];
%         n = n+1;
% %         end
%     end
% 
% end
% 
% 
% figure(55);
% scatter3(data(:,1),data(:,2),data(:,3));
% title('3D Point Cloud');
% xlabel('X (m)');
% ylabel('Y (m)');
% 
% 
% %  step = 0.05;
% 
% xgrid = min(x):step:max(x)+1;
% ygrid = min(y):step:max(y)+1;
% zgrid = min(z):step:max(z)+1;
% 
% xsize  = size(xgrid,2)+1;
% ysize  = size(ygrid,2)+1;
% zsize  = size(zgrid,2)+1;
% 
% xori = min(x);
% yori = min(y);
% zori = min(z);
% 
% Z = zeros(xsize,ysize);
% 
% for k=1:length(x)
% %     
% %     kx = round(Data(k,1)*100) + xsize;
% %     ky = round(Data(k,2)*100) + ysize;
% 
%     kx = round((x(k,1)-xori)/step)+1;
%     ky = round((y(k,1)-yori)/step)+1;
%     Z(kx,ky) = Z(kx,ky) +1;
%     
% end
% 
% figure;
% imshow(Z);
% 




% load seamount 
% dat = [-y,x]; % Grid corrected for negative y-values
% hold on 
% hist3(dat) % Draw histogram in 2D 
% 
% n = hist3(dat); % Extract histogram data;
%                 % default to 10x10 bins
% n1 = n'; 
% n1( size(n,1) + 1 ,size(n,2) + 1 ) = 0; 
% 
% 
% xb = linspace(min(dat(:,1)),max(dat(:,1)),size(n,1)+1);
% yb = linspace(min(dat(:,2)),max(dat(:,2)),size(n,1)+1);
% 
% 
% h = pcolor(xb,yb,n1);
% 
% 
% set(h, 'zdata', ones(size(n1)) * -max(max(n))) 
% colormap(hot) % heat map 
% title('Seamount:Data Point Density Histogram and Intensity Map');
% grid on 
% 
% view(3);
% 
% load carbig
% X = [MPG,Weight];
% hist3(X,[7 7]);
% xlabel('MPG'); ylabel('Weight');
% set(gcf,'renderer','opengl');
% set(get(gca,'child'),'FaceColor','interp','CDataMode',...
% 'auto');


