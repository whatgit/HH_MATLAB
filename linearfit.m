function [rl ll] = linearfit(cloudpoint)

% 
% cloud = load('X:\MT\Mapfilterstuff\output\pointcloud2.txt'); 
% 
% step = 0.05
% [data r c] = transOrth(cloud,step); % Rotate Data till Dominant Direction Arg: points,step
% 
% data = OutlierRej(data); % Reject Outliers
% 
% [cloudpoint r c] = transOrth(data,step); % Return Data Rotated along dominant Direction. Arg: points,step

% Seperation
lidx = 1;
ridx = 1;
left = [];
right = [];

for kk = 1:length(cloudpoint)
    if(cloudpoint(kk,1)<0)
        left(lidx,:) = [cloudpoint(kk,1) cloudpoint(kk,2) cloudpoint(kk,3)];
        lidx = lidx + 1;
    else
        right(ridx,:) = [cloudpoint(kk,1) cloudpoint(kk,2) cloudpoint(kk,3)];
        ridx = ridx + 1;
    end

end

pl = polyfit(left(:,1),left(:,2),1)
xsaml = -3.5:0.1:-0;
fl = polyval(pl,xsaml);

figure;
plot(xsaml,fl,'-');hold on;
axis equal


pr = polyfit(right(:,1),right(:,2),1)
xsamr = -5:0.1:15;
fr = polyval(pr,xsamr);

plot(xsamr,fr,'-');hold on;
axis equal


scatter(left(:,1),left(:,2),'r','.');hold on;
scatter(right(:,1),right(:,2),'r','.');hold on;
title('2D Point Cloud');
xlabel('X (m)');
ylabel('Y (m)');
axis equal


Fleft = mean(left(:,1));
Fright = mean(right(:,1));
leftline = min(left(:,2)):0.1: max(left(:,2));
rightline = min(right(:,2)):0.1: max(right(:,2));

for kk = 1:length(leftline)
    ll(kk,:) = [Fleft leftline(kk)];
end
    
for kk = 1:length(rightline)
    rl(kk,:) = [Fright rightline(kk)];
end


figure;
scatter(left(:,1),left(:,2),'r','.');hold on;
scatter(right(:,1),right(:,2),'r','.');hold on;

plot(rl(:,1),rl(:,2),'-','LineWidth',4);hold on;
plot(ll(:,1),ll(:,2),'-','LineWidth',4);hold on;

title('Line Fitting');
xlabel('X (m)');
ylabel('Y (m)');
axis equal





end