function [pillarcloud pillarL pillarR MeanDleft MeanDright] = refinePillar(cloudpoint,leftp,rightp)

% Function   : refinePillar
% 
% Purpose    : Combine Pillar estmation and pointcould data to get a joint
%              estimation
% 
% Parameters : cloudpoint, Point cloud data
%              leftp and rightp, estimated pillar position
% 
% Return     : pillarcloud, candidates that belongs to a pillar in point
%              cloud
%                pillarL, Joint estimated pillar on left side
%                pillarR, Joint estimated pillar on right side
%                MeanDleft, Mean value of left pillar
%                MeanDright, Mean value of right pillar



n=1;
Data = [];

x = cloudpoint(:,1);
y = cloudpoint(:,2);
z = cloudpoint(:,3);

 %%
 
 leftp;
 rightp;
 x;
 y;
 z;
 
meanl= [];
varl= [];
meanr= [];
varr= [];
 
pcleft = [];
leftidx = 1;
pcright = [];
rightidx = 1;
 
rpcloud =[];
rpcidx = 1;
sidx = 1;

 tmp = [0 0 0];
  
 for k = 1:length(leftp)
     curp = leftp(k,:);
     for kk = 1:length(x)
%         dist = sqrt((curp(1,1)-x(kk))^2+(curp(1,2)-y(kk))^2);
        if abs(curp(1,2)-y(kk))<0.55 && (x(kk)<=0) && ((x(kk)>-5))
            rpcloud(rpcidx,:) = [x(kk) y(kk) z(kk) k];
            rpcidx = rpcidx +1;
            tmp(sidx,:) =[x(kk) y(kk) z(kk)];
            sidx = sidx +1;
            % Store to left pillar
            pcleft(leftidx,:) = [x(kk) y(kk) z(kk) k];
            leftidx = leftidx +1 ;
        end
        
        if kk == length(x)
            meanl(k,:) =  [mean(tmp(:,1)) mean(tmp(:,2)) mean(tmp(:,3))]
            varl(k,:) = [var(tmp(:,1)) var(tmp(:,2)) var(tmp(:,3))];
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
%       if abs(curp(1,2)-y(kk))<0.55 && (x(kk)>=0) && ((x(kk)<5))
           if abs(curp(1,2)-y(kk))<0.55 && (x(kk)>=0) && ((x(kk)<5))
            rpcloud(rpcidx,:) = [x(kk) y(kk) z(kk) k];
            rpcidx = rpcidx +1;
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

  pillarcloud = rpcloud;

  pillarL = meanl;
  pillarR = meanr;
  pillarL(:,3) = zeros(length(meanl),1);
  pillarR(:,3) = zeros(length(meanr),1);
  
%%
% figure(50);
% scatter3(rpcloud(:,1), rpcloud(:,2), rpcloud(:,3),'.');


lx1 = zeros(size(leftp(:,1)))
lx1 = lx1-3
lx2 = lx1+2.8

rx1 = zeros(size(rightp(:,1)));
rx1 = rx1+3;
rx2 = rx1-2.8;

ly1 = zeros(size(leftp(:,2)));
ly1 = leftp(:,2)-0.55;
ly2 = ly1+1.1;


ry1 = zeros(size(rightp(:,2)));
ry1 = rightp(:,2)-0.55;
ry2 = ry1+1.1;


figure(51111);
% scatter(rpcloud(:,1), rpcloud(:,2)); hold on; 
scatter(cloudpoint(:,1), cloudpoint(:,2)); hold on; 
h=plot(meanl(:,1),meanl(:,2),'.'); hold on;
set(h,'MarkerSize',25,'color','red');
h=plot(meanr(:,1),meanr(:,2),'.');hold on;  
set(h,'MarkerSize',25,'color','red');
for nn = 1:length(lx1)
plot([lx1(nn),lx2(nn)],[ly1(nn),ly1(nn)],'r');hold on;
plot([lx1(nn),lx2(nn)],[ly2(nn),ly2(nn)],'r');hold on;
plot([lx1(nn),lx1(nn)],[ly1(nn),ly2(nn)],'r');hold on
plot([lx2(nn),lx2(nn)],[ly1(nn),ly2(nn)],'r');hold on;

end

for nn = 1:length(rx1)
plot([rx1(nn),rx2(nn)],[ry1(nn),ry1(nn)],'r');hold on;
plot([rx1(nn),rx2(nn)],[ry2(nn),ry2(nn)],'r');hold on;
plot([rx1(nn),rx1(nn)],[ry1(nn),ry2(nn)],'r');hold on
plot([rx2(nn),rx2(nn)],[ry1(nn),ry2(nn)],'r');hold on;

end

title('Pillar Estimated');
xlabel('X (m)');
ylabel('Y (m)');
axis equal

hold off;



figure(33333);
scatter(cloudpoint(:,1), cloudpoint(:,2)); hold on; 
h=plot(leftp(:,1),leftp(:,2),'.'); hold on;
set(h,'MarkerSize',20,'color','green');
h=plot(rightp(:,1),rightp(:,2),'.');hold on;  
set(h,'MarkerSize',20,'color','green');
h=plot(meanl(:,1),meanl(:,2),'.'); hold on;
set(h,'MarkerSize',15,'color','red');
h=plot(meanr(:,1),meanr(:,2),'.');hold on;  
set(h,'MarkerSize',15,'color','red');

for nn = 1:length(lx1)
plot([lx1(nn),lx2(nn)],[ly1(nn),ly1(nn)],'r');hold on;
plot([lx1(nn),lx2(nn)],[ly2(nn),ly2(nn)],'r');hold on;
plot([lx1(nn),lx1(nn)],[ly1(nn),ly2(nn)],'r');hold on
plot([lx2(nn),lx2(nn)],[ly1(nn),ly2(nn)],'r');hold on;

end

for nn = 1:length(rx1)
plot([rx1(nn),rx2(nn)],[ry1(nn),ry1(nn)],'r');hold on;
plot([rx1(nn),rx2(nn)],[ry2(nn),ry2(nn)],'r');hold on;
plot([rx1(nn),rx1(nn)],[ry1(nn),ry2(nn)],'r');hold on
plot([rx2(nn),rx2(nn)],[ry1(nn),ry2(nn)],'r');hold on;

end

title('Pillar Estimated');
xlabel('X (m)');
ylabel('Y (m)');
axis equal

hold off;

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
MeanDleft = [mean(leftd(:,1)) mean(leftd(:,2))]

varl = var(leftd)
for r = 1:length(meanr)-1
    rightd(r,:) = [meanr(r+1,1)-meanr(r,1) meanr(r+1,2)-meanr(r,2)];
    
end

MaxDright = [max(rightd(:,1)) max(rightd(:,2))]
MeanDright= [mean(rightd(:,1)) mean(rightd(:,2))]
varr = var(rightd)

end