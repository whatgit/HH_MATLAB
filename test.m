% This program generate a joint estimation of a pillar map
% Input: 1.Point cloud data
%        2.Pillar position data
% Output: 1.Pillar map in 2D metric scale
%         2.Mean distance between pillars and corridor width
%

close all;clear all;

% Data file, modify the path below
% clouddata <= point cloud data
% pillarPosdata <= Pillar position data
clouddata = load('X:\MT\Mapfilterstuff\output\pointcloud2.txt'); 
pillarPosdata = load('X:\MT\Mapfilterstuff\output\warehouse4_houghline2.txt');


% Data preprocessing
%
step = 0.05; % Step of mesh grid for Randon Transform
data = transOrth(clouddata,step); % Find dominant orientation, data: rotated point cloud
data = OutlierRej(data); % Reject Outliers and tags
cloudpoint = transOrthxy(data,step); % Find dominant orientation without outliers and tags
[leftp rightp] = estimatepillar(pillarPosdata); % Return estimated Pillar according to ImgProcess

% Manually Cut off some Pillars
for k = 1:(length(leftp)-2)
    cutleftp(k,:)=leftp(k,:); 
end

for k = 1:(length(rightp)-3)
    cutrightp(k,:)=rightp(k,:); 
end
 
figure(1);
scatter(cloudpoint(:,1),cloudpoint(:,2));
title('2D Point Cloud');
xlabel('X (m)');
ylabel('Y (m)');
axis equal
 
% Data Fusion
[pillarcloud pillarL pillarR MeanDleft MeanDright] = refinePillar(cloudpoint,cutleftp,cutrightp); % Refine according to pointcloud
[Distance mid xleft xright]= findwidth(pillarcloud); % Find width of corridor
[pl pr Corrwidth] = pillarfromimgp(pillarL,pillarR,xright,xleft,pillarcloud); % replace x position with corridor width
[MeanDleft MeanDright varl varr]=calmeanvar(pillarL,pillarR);% Calculate means and vars


% [rl ll] = linearfit(cloudpoint);
% [pl pr Corrwidth] = pillarfromimgp(leftp,rightp,xright,xleft,pillarcloud);
% [pl pr Corrwidth] =-+ pillarfromimgp(pillarL,pillarR,rl,ll,pillarcloud);

figure(1);
scatter(cloudpoint(:,1),cloudpoint(:,2));
title('2D Point Cloud');
xlabel('X (m)');
ylabel('Y (m)');
axis equal


 

%%

figure;
scatter(data(:,1),data(:,3));
title('2D Point Cloud');
xlabel('X (m)');
ylabel('Y (m)');
axis equal

figure;
scatter3(data(:,1),data(:,2),data(:,3));
title('3D Point Cloud');
xlabel('X (m)');
ylabel('Y (m)');
axis equal

figure;
scatter3(pillarcloud(:,1),pillarcloud(:,2),pillarcloud(:,3),'.');hold on;
scatter3(pillarL(:,1),pillarL(:,2),pillarL(:,3),'o','r','filled');hold on;
scatter3(pillarR(:,1),pillarR(:,2),pillarR(:,3),'o','r','filled');hold on;
title('3D Pillar Cloud');
xlabel('X (m)');
ylabel('Y (m)');
axis equal

% NA = plotPillar(pillarL,pillarR,pillarcloud);
 NA = plotPillar(pl,pr,pillarcloud);
 
 scaleLength = (abs(MeanDleft(1,2)) + abs(MeanDright(1,2)))*0.5/2.70
 scaleWidth = abs(Distance)/2.95