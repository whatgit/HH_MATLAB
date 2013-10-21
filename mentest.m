
cloud = load('X:\MT\Mapfilterstuff\output\pointcloud2.txt'); 

step = 0.05;
[data r c] = transOrth(cloud,step); % Rotate Data till Dominant Direction Arg: points,step

xx = data(:,1);
yy = data(:,2);
zz = data(:,3);
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

step = 0.01

xgrid = min(x):step:max(x)+1;
ygrid = min(y):step:max(y)+1;
% zgrid = min(z):step:max(z)+1;

xsize  = size(xgrid,2)+1;
ysize  = size(ygrid,2)+1;
% zsize  = size(zgrid,2)+1;

xori = min(x);
yori = min(y);
% zori = min(z);
Z = zeros(xsize,ysize);
for k=1:length(x)
%     
%     kx = round(Data(k,1)*100) + xsize;
%     ky = round(Data(k,2)*100) + ysize;

    kx = round((x(k,1)-xori)/step)+1;
    ky = round((y(k,1)-yori)/step)+1;
    Z(kx,ky) = Z(kx,ky) +1;
    
end

% imshow(Z);

% CorrWidth = 2.90;
% PillarWidth = 2.70;
% step = 0.01;
% % Pnum = length(rightp,1);
% Pnum = 12;
% 
% xgrid = -0.2-CorrWidth/2:step:CorrWidth+0.2;
% ygrid = -0.2-PillarWidth*Pnum:step:PillarWidth*Pnum+0.2;



% mu = [round(CorrWidth/(step*2)) round(CorrWidth/(step*2))];
% Sigma = [.25 .3; .3 1];

mu = [0 0]

xmesh = -0.15:step:0.15;
ymesh = -0.15:step:0.15;

[X1,X2] = meshgrid(xmesh,ymesh);
F = mvnpdf([X1(:) X2(:)],mu);
F = reshape(F,length(ymesh),length(xmesh));
figure(3);
surf(xmesh,ymesh,F);
caxis([min(F(:))-.5*range(F(:)),max(F(:))]);
% axis([-3 3 -3 3 0 .4])
axis equal
xlabel('x1'); ylabel('x2'); zlabel('Probability Density');







