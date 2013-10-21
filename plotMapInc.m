% function plotMapInc
% tagdata = load('D:\Warehouse\Pilarpose1.txt');
% tagdata = load('X:\MT\MapVisual\tagnew.txt')
% tagdata = load('D:\Thesis\Warehouse Bag\Lab\output\gym2_3_KF_map_Thresh100px.txt');
clear all;
close all;
% tagdata = load('/media/fyt/Os2/MT/Mapfilterstuff/output/warehouse4_filter6.txt');

tagdata = load('X:\MT\Mapfilterstuff\output\warehouse4_pointcloud4.txt');

pos = [tagdata(:,2) tagdata(:,3)];
z = tagdata(:,4);
heading = tagdata(:,5);

%% Set Para
% Mcx = 0.871576; 
% Mcy = 0.493151;

Mcx = 16/sqrt(337);
Mcy = 9/sqrt(337);

% alfa = pi/2;
% R = [cos(alfa) -sin(alfa); sin(alfa) cos(alfa)];
% 
% posNew = pos * R;
% tagNew = tag * R;

%% Cal
ppxl1 = [tagdata(:,7) tagdata(:,8)];
ppxl2 = [tagdata(:,9) tagdata(:,10)];

% pwc1(:,1) =((tagdata(:,5)-320)/640).*(Mcx.*z.*2);
% pwc1(:,2) =(-(tagdata(:,6)-160)/320).*(Mcy.*z.*2);
% pwc2(:,1) =((tagdata(:,7)-320)/640).*(Mcx.*z.*2);
% pwc2(:,2) =(-(tagdata(:,8)-160)/320).*(Mcy.*z.*2);

    pwc1(:,1) = ((tagdata(:,7)-320.0).*((2.0*tagdata(:,4))/734.3));
    pwc1(:,2) = ((180.0-tagdata(:,8)).*((2.0*tagdata(:,4))/734.3));
    pwc2(:,1) = ((tagdata(:,9)-320.0).*((2.0*tagdata(:,4))/734.3));
    pwc2(:,2) = ((180-tagdata(:,10)).*((2.0*tagdata(:,4))/734.3));

% pwc1(:,1) =tagdata(:,2)+((tagdata(:,5)-320)/640).*(Mcx.*z.*2);
% pwc1(:,2) =tagdata(:,3)+((tagdata(:,6)-160)/320).*(Mcy.*z.*2);
% pwc2(:,1) =tagdata(:,2)+((tagdata(:,7)-320)/640).*(Mcx.*z.*2);
% pwc2(:,2) =tagdata(:,3)+((tagdata(:,8)-160)/320).*(Mcy.*z.*2);

% pos1 = [pwc1(:,1) pwc1(:,2)] * R;
% pos2 = [pwc2(:,1) pwc2(:,2)] * R;

%% Rotate 

iIndex = size(heading,1);

for i=1:iIndex
    alfa = heading(i,1)*3.14/180;
    R = [cos(alfa) -sin(alfa); sin(alfa) cos(alfa)];
%     pos1(i,:) = pwc1(i,:) * R;
%     pos2(i,:) = pwc2(i,:) * R;
    pos1(i,:) = pwc1(i,:) ;
    pos2(i,:) = pwc2(i,:) ;
    
end

%% Convert to WF
pwc1(:,1) =tagdata(:,2)+pos1(:,1);
pwc1(:,2) =tagdata(:,3)+pos1(:,2);
pwc2(:,1) =tagdata(:,2)+pos2(:,1);
pwc2(:,2) =tagdata(:,3)+pos2(:,2);
    
    





%% Plot Pillar 
    figure;
    
    lidx = 1;
    ridx = 1;
    for i = 1:iIndex

%     plot(pwc1(i,1),pwc1(i,2),'--rs','LineWidth',0.1);hold on;
%                    
%                    plot(pwc2(i,1),pwc2(i,2),'--rs','LineWidth',0.1); hold on;
%     
    plot([pwc1(i,1) pwc2(i,1)],[pwc1(i,2) pwc2(i,2)]);hold on;

    if( (pwc1(i,1) < pos(i,1)) && (pwc2(i,1) < pos(i,1)) )
       left(lidx,:) =  [pwc2(i,1) pwc2(i,2)];
       lidx = lidx +1;
    end
    
    if( (pwc1(i,1) > pos(i,1)) && (pwc2(i,1) > pos(i,1)) )
       right(ridx,:) =  [pwc1(i,1) pwc1(i,2)];
       ridx = ridx +1;
    end
        
    end
    
scatter(pos(:,1),pos(:,2),'r');
axis equal
title('Map');
xlabel('X (m)');
ylabel('Y (m)');
    
%% Plot Trace
    
figure;
  scatter(pwc1(:,1),pwc1(:,2),'b'); hold on;
  scatter(pwc2(:,1),pwc2(:,2),'g');hold on;
%  scatter(pwc1(:,1),pwc1(:,2),'b'); hold on;
%  scatter(pwc2(:,1),pwc2(:,2),'g');hold on;
 scatter(pos(:,1),pos(:,2),'r');
 axis equal
 
% scatter(tag(:,1),tag(:,2));
% end

figure;
  scatter(left(:,1),left(:,2),'g'); hold on;
  scatter(right(:,1),right(:,2),'b');hold on;
%  scatter(pwc1(:,1),pwc1(:,2),'b'); hold on;
%  scatter(pwc2(:,1),pwc2(:,2),'g');hold on;
 scatter(pos(:,1),pos(:,2),'r');
 axis equal
 
idxss11 = 1;
for idxss = 1: length(right)
    if(abs(right(idxss,1))>0.3)
        abc(idxss11,:) = right(idxss,:)
        idxss11 = idxss11+1;
    end
    
end

right = [];
right = abc;
 
 %%
threshold = 1.2;
tmpx = 0;
tmpy = 0;
cnt = 1;
lidx = 1;
ridx = 1;
leftpillar = [];
rightpillar = [];


leftwidx = [];
nbr = 1;
leftwidx(1) = nbr

for i = 1:(length(left)-1)
 if( abs(left(i,2) - left(i+1,2)) >= threshold )
     nbr = nbr + 1;
     leftwidx(i+1) = nbr;
 end
 
 if( abs(left(i,2) - left(i+1,2)) < threshold )
     leftwidx(i+1) = nbr;
 end
 
end

rightwidx = [];
 nbr = 1;
rightwidx(1) = nbr
for i = 1:(length(right)-1)
 if( abs(right(i,2) - right(i+1,2)) >= threshold )
     nbr = nbr + 1;
     rightwidx(i+1) = nbr;
 end
 
 if( abs(right(i,2) - right(i+1,2)) < threshold )
     rightwidx(i+1) = nbr;
 end
 
end
 
n = 1; 

 for i = 1:(length(left))
     
     if leftwidx(i)== n
         tmpy = left(i,2) + tmpy;
         tmpx = left(i,1) + tmpx;
         cnt = cnt+1; 
     end
     
     if leftwidx(i) ~= n
        leftp(n,:) = [tmpx/cnt tmpy/cnt];
        cnt = 1;
        n = n+1;
         tmpy = left(i,2);
         tmpx = left(i,1);      
     end
    
     if i == (length(left))
         if leftwidx(i) == leftwidx(i-1)
             leftp(n,:) = [tmpx/cnt tmpy/cnt];         
         end  
     end
 end
 
 tmpx = 0;
tmpy = 0;
cnt = 1;
 
 
n = 1; 

 for i = 1:(length(right))
     
     if rightwidx(i)== n
         tmpy = right(i,2) + tmpy;
         tmpx = right(i,1) + tmpx;
         cnt = cnt+1; 
     end
     
     if rightwidx(i) ~= n
        rightp(n,:) = [tmpx/cnt tmpy/cnt];
        cnt = 1;
        n = n+1;
         tmpy = right(i,2);
         tmpx = right(i,1);      
     end
    
     if i == (length(right))
         if rightwidx(i) == rightwidx(i-1)
             rightp(n,:) = [tmpx/cnt tmpy/cnt];         
         end  
     end
 end
 
 
%  figure(10);
% %   scatter(leftp(:,1),leftp(:,2),'g'); hold on;
% %   scatter(rightp(:,1),rightp(:,2),'b');hold on;
%   h=plot(leftp(:,1),leftp(:,2),'.'); hold on;
%   set(h,'MarkerSize',25,'color','red');
%   h=plot(rightp(:,1),rightp(:,2),'.');hold on;  
% set(h,'MarkerSize',25,'color','red');
% %  scatter(pwc1(:,1),pwc1(:,2),'b'); hold on;
% %  scatter(pwc2(:,1),pwc2(:,2),'g');hold on;
%  scatter(pos(:,1),pos(:,2),'r'); hold on;
%  axis([-20 20 -5 35])
 
 %% Analysis
 
%  leftlength = length(leftp);
%  rightlength = length(rightp);
%  
%  leftAvr = sum(leftp(:,1))/length(leftp);
%  
%  
%  p = polyfit(leftp(:,1),leftp(:,2),1)
%  
%   x = leftp(:,1)
% %   x= linspace(0,leftp(length(leftp),1), n); % Adapt n for resolution of graph
%   y= p(1,1)+p(1,2)*x;
% %   figure(5),
%   plot(leftp(:,1),leftp(:,2),'o'); hold on;
%   plot(x,y);
%  
%  
%%
 tagdata = load('X:\MT\Mapfilterstuff\output\pointcloud2.txt');

n=1;
Data = [];


thetayz = -4.8;
Ryz = [cosd(thetayz) -sind(thetayz); sind(thetayz) cosd(thetayz)];

 for k=1:length(tagdata(:,3))
%     tagdata(k,3) = tagdata(k,3) + 0.42;
    
    vv = [tagdata(k,3) tagdata(k,2)] * Ryz;
    tagdata(k,3) = vv(1,1);
    tagdata(k,2) = vv(1,2);
%         plot3(x(k), y(k), z(k), 'ws--', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
 end
 
thetaxy = 2.5;
Rxy = [cosd(thetaxy) -sind(thetaxy); sind(thetaxy) cosd(thetaxy)];
 
  for k=1:length(tagdata(:,1))
%     tagdata(k,3) = tagdata(k,3) + 0.42;
    
    vv = [tagdata(k,2) tagdata(k,1)] * Rxy;
    tagdata(k,2) = vv(1,1);
    tagdata(k,1) = vv(1,2);
%         plot3(x(k), y(k), z(k), 'ws--', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
 end


 x = tagdata(:,1);
y = tagdata(:,2);
z = tagdata(:,3);

% Get rid of noise

for kk = 1:length(x)
    if((abs(x(kk))<5)&&((abs(y(kk))<100))&&((z(kk)<5.5)&&(z(kk)>-1.2)))
        if ((x(kk)<-0.4)||(abs(y(kk))>23)||(x(kk))>1)
        Data(n,1) = x(kk);
        Data(n,2) = y(kk);
        Data(n,3) = z(kk);
        n = n+1;
        end
    end

end

x = Data(:,1);
y = Data(:,2);
z = Data(:,3);

 %%
 
 leftp;
 rightp;
 x;
 y;
 z;
 pos;
 
pcleft = [];
leftidx = 1;
pcright = [];
rightidx = 1;
meanl= [];
varl= [];
meanr= [];
varr= [];
 
rpcloud =[];
 rpcidx = 1;
 sidx = 1;

 tmp = [0 0 0];
  
 for k = 1:length(leftp)
     curp = leftp(k,:);
     for kk = 1:length(x)
%         dist = sqrt((curp(1,1)-x(kk))^2+(curp(1,2)-y(kk))^2);
        if abs(curp(1,2)-y(kk))<0.55 && (x(kk)<=1) && ((x(kk)>-4))
            rpcloud(rpcidx,:) = [x(kk) y(kk) z(kk) k];
            rpcidx = rpcidx +1
            tmp(sidx,:) =[x(kk) y(kk) z(kk)];
            sidx = sidx +1;
            % Store to left pillar
            pcleft(leftidx,:) = [x(kk) y(kk) z(kk) k];
            leftidx = leftidx +1 ;
        end
        
        if kk == length(x)
            meanl(k,:) =  [mean(tmp(:,1)) mean(tmp(:,2)) mean(tmp(:,3))]
            varl(k,:) = [var(tmp(:,1)) var(tmp(:,2)) var(tmp(:,3))]
             sidx = 1;
              tmp = [0 0 0];
        end
         
     end
          
 end
 
 
  for k = 1:length(rightp)
     curp = rightp(k,:);
     
     for kk = 1:length(x)
%         dist = sqrt((curp(1,1)-x(kk))^2+(curp(1,2)-y(kk))^2);
%         if dist < 0.8
      if abs(curp(1,2)-y(kk))<0.55 && (x(kk)>=1) && ((x(kk)<8))
            rpcloud(rpcidx,:) = [x(kk) y(kk) z(kk) k];
            rpcidx = rpcidx +1
            tmp(sidx,:) =[x(kk) y(kk) z(kk)];
            sidx = sidx +1;
            % Store to right pillar
            pcright(rightidx,:) = [x(kk) y(kk) z(kk) k];
            rightidx = rightidx +1 ;
            
      end
        
       if kk == length(x)
            meanr(k,:) =  [mean(tmp(:,1)) mean(tmp(:,2)) mean(tmp(:,3))]
            varr(k,:) = [var(tmp(:,1)) var(tmp(:,2)) var(tmp(:,3))]
             sidx = 1;
            tmp = [0 0 0];
        end
         
     end
          
  end


 %%
 
%  theta = -4;
% R = [cosd(theta) -sind(theta); sind(theta) cosd(theta)];
%  
%  
%  for k=1:length(rpcloud)
%     rpcloud(k,3) = rpcloud(k,3) + 0.42;
%     
%     vv = [rpcloud(k,3) rpcloud(k,2)] * R;
%     rpcloud(k,3) = vv(1,1);
%     rpcloud(k,2) = vv(1,2);
% %         plot3(x(k), y(k), z(k), 'ws--', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
%  end
% 
%   for k=1:length(meanl)
%     meanl(1,3) = meanl(1,3) + 0.42;
%     
%     vv = [meanl(k,3) meanl(k,2)] * R;
%     meanl(k,3) = vv(1,1);
%     meanl(k,2) = vv(1,2);
% %         plot3(x(k), y(k), z(k), 'ws--', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
%   end
% 
%  for k=1:length(meanr)
%     meanr(1,3) = meanr(1,3) + 0.42;
%     
%     vv = [meanr(k,3) meanr(k,2)] * R;
%     meanr(k,3) = vv(1,1);
%     meanr(k,2) = vv(1,2);
% %         plot3(x(k), y(k), z(k), 'ws--', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
% end
%  
 
 

%%
figure(50);
scatter3(rpcloud(:,1), rpcloud(:,2), rpcloud(:,3),'.');

figure(51);
scatter(rpcloud(:,1), rpcloud(:,2)); hold on; 
  h=plot(meanl(:,1),meanl(:,2),'.'); hold on;
  set(h,'MarkerSize',25,'color','red');
  h=plot(meanr(:,1),meanr(:,2),'.');hold on;  
set(h,'MarkerSize',25,'color','red');


hold off;
% data = vdhf_data.data(:,4);
% scatter3(x, y, z);
% [xx,yy] = meshgrid(x, y);
% scatter(x,y)

%%

% figure(50);
% for k=1:length(rpcloud)
%     scatter3(rpcloud(k,1), rpcloud(k,2), rpcloud(k,3),'.'); hold on;    
%    k
%     
% %         plot3(x(k), y(k), z(k), 'ws--', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
% end
% 
% 
% figure(51);
% for k=1:length(rpcloud)
%     scatter(rpcloud(k,1), rpcloud(k,2)); hold on;  
% 
%    k
%     
% %         plot3(x(k), y(k), z(k), 'ws--', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
% end
% 
%   h=plot(meanl(:,1),meanl(:,2),'.'); hold on;
%   set(h,'MarkerSize',25,'color','red');
%   h=plot(meanr(:,1),meanr(:,2),'.');hold on;  
% set(h,'MarkerSize',25,'color','red');
% 
% 
% hold off;





% 
% 
% %%
% figure(10);hold on;
% for k=1:length(Data)
%     scatter(x(k), y(k),'.'); hold on;
%     k
% %         plot3(x(k), y(k), z(k), 'ws--', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
% end
% title('Pillar Point Cloud');
% xlabel('X');
% ylabel('Y');
% 
%  figure(10);
% %   scatter(leftp(:,1),leftp(:,2),'g'); hold on;
% %   scatter(rightp(:,1),rightp(:,2),'b');hold on;
%   h=plot(leftp(:,1),leftp(:,2),'.'); hold on;
%   set(h,'MarkerSize',25,'color','green');
%   h=plot(rightp(:,1),rightp(:,2),'.');hold on;  
% set(h,'MarkerSize',25,'color','green');
% %  scatter(pwc1(:,1),pwc1(:,2),'b'); hold on;
% %  scatter(pwc2(:,1),pwc2(:,2),'g');hold on;
%  scatter(pos(:,1),pos(:,2),'r'); hold on;
%  axis([-20 20 -5 35]);
%  hold off;
% 
% 
% 
% 
% %%
% figure(12);hold on;
% for k=1:length(Data)
%     scatter(z(k), y(k)); hold on;
%    k
%     
% %         plot3(x(k), y(k), z(k), 'ws--', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
% end
% 
%       kk=plot(leftp(:,1)-leftp(:,1),leftp(:,2),'.'); hold on;
%       set(kk,'MarkerSize',25,'color','red');
%       kk=plot(rightp(:,1)-rightp(:,1),rightp(:,2),'.');hold on;  
%       set(kk,'MarkerSize',25,'color','green');
% 
% title('Pillar Point Cloud');
% xlabel('Z');
% ylabel('Y');
% 
%  hold off;
%  
%  
% %%
% figure(50);
% for k=1:length(rpcloud)
%     scatter3(rpcloud(k,1), rpcloud(k,2), rpcloud(k,3)); hold on;    
%    k
%     
% %         plot3(x(k), y(k), z(k), 'ws--', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
% end
% 
% %%
% 
% figure(51);
% for k=1:length(rpcloud)
%     scatter(rpcloud(k,1), rpcloud(k,2)); hold on;    
%    k
%     
% %         plot3(x(k), y(k), z(k), 'ws--', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
% end
% 
% %%
% figure(53);
% for k=1:length(rpcloud)
%     scatter(rpcloud(k,3), rpcloud(k,2)); hold on;    
%    k
%     
% %         plot3(x(k), y(k), z(k), 'ws--', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
% end
% 
% %%
% 


%% Calculate d & m



pcleft;
pcright;
meanl;
meanr;
leftd = [];
rightd = [];

for l = 1:length(meanl)-1
    leftd(l,:) = [meanl(l+1,1)-meanl(l,1) meanl(l+1,2)-meanl(l,2)];
    
end

MaxDleft = [max(leftd(:,1)) max(leftd(:,2))]
MeanDeft = [mean(leftd(:,1)) mean(leftd(:,2))]

varll = var(leftd);
for r = 1:length(meanr)-1
    rightd(r,:) = [meanr(r+1,1)-meanr(r,1) meanr(r+1,2)-meanr(r,2)];
    
end

MaxDright = [max(rightd(:,1)) max(rightd(:,2))]
MeanDright= [mean(rightd(:,1)) mean(rightd(:,2))]
varrr = var(rightd);

%%
% xgird = -5:0.01:4.99;
% ygird = -3:0.01:32.99;
% xori = -5;
% yori = -3;
% xgird = -5:0.01:4.99;
% ygird = -3:0.01:32.99;
% 
% [xx,yy]= meshgrid(xgird,ygird);

% Z = zeros(1000,3600);
% ZZ = zeros(1001,3601);
% step = 0.01;


xgrid = -5:0.05:4.99;
ygrid = -3:0.05:32.99;
xsize  = size(xgrid,2);
ysize  = size(ygrid,2);

xori = -5;
yori = -3;

Z = zeros(xsize,ysize);

for k=1:length(rpcloud(:,1))
%     
%     kx = round(Data(k,1)*100) + xsize;
%     ky = round(Data(k,2)*100) + ysize;

    kx = round((rpcloud(k,1)-xori)/0.05);
    ky = round((rpcloud(k,2)-yori)/0.05);
    Z(kx,ky) = Z(kx,ky) +1;
    
end

theta = 0:0.1:180;
[R,xp] = radon(Z,theta);
imshow(R,[],'Xdata',theta,'Ydata',xp,...
            'InitialMagnification','fit')
xlabel('\theta (degrees)')
ylabel('x''')
colormap(hot), colorbar
iptsetpref('ImshowAxesVisible','off')

[r,c] = find (R == max(max(R)))

imshow(imrotate(R,(180-c)))


