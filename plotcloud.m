clc; clear all; close all


%%
% tagdata = load('/media/fyt/Os2/MT/Mapfilterstuff/output/pointcloud2.txt');
tagdata = load('X:\MT\Mapfilterstuff\output\pointcloud2.txt');

x = tagdata(:,1);
y = tagdata(:,2);
z = tagdata(:,3);
n=1;
Data = [];

% Get rid of noise

for kk = 1:length(x)
    if((abs(x(kk))<5)&&((abs(y(kk))<100))&&((z(kk)<5.5)&&(z(kk)>-1.2)))
        Data(n,1) = x(kk);
        Data(n,2) = y(kk);
        Data(n,3) = z(kk);
        n = n+1;
    end

end


x = Data(:,1);
y = Data(:,2);
z = Data(:,3);
% data = vdhf_data.data(:,4);
% scatter3(x, y, z);
% [xx,yy] = meshgrid(x, y);
% scatter(x,y)

%%
figure(1);hold on;
for k=1:length(Data)
    scatter(x(k), y(k)); hold on;
    
%         plot3(x(k), y(k), z(k), 'ws--', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
end
title('Pillar Point Cloud');
xlabel('X');
ylabel('Y');
hold off;

%%
figure(2);hold on;
for k=1:length(Data)
    scatter(z(k), y(k)); hold on;
    
%         plot3(x(k), y(k), z(k), 'ws--', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
end
title('Pillar Point Cloud');
xlabel('Z');
ylabel('Y');
hold off;


%%


%%
% figure(10);hold on;
% for k=1:2000
%     scatter(z(k), y(k)); hold on;
%     
% %         plot3(x(k), y(k), z(k), 'ws--', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
% 
% end
% 
% 
% 
% %%
% 
% 
% figure(10);hold on;
% for k=1:6000
%     scatter(z(k), y(k)); hold on;
%     
% %         plot3(x(k), y(k), z(k), 'ws--', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
% 
% end
% 
% 
% 
% 
% 
% %% generates the first figure using 'plot3'
% figure
% plot3(x,y,z)
% grid on
% 
% %% generates the second figure using 'meshc' to include the
% % contour in the figure, and rotates the figure with 'view'
% figure
% meshc(x,y,z)
% view(-37, 15)
% 
% %% generates the third 3D figure using 'surfc' to include the
% % contour in the image, and also rotates the figure with 'view'
% figure
% surfc(x,y,z)
% view(-47, 25)
% 
% 
% % 
% % grid on
% % hold all
% % 
% % for k=1:length(x)
% % 
% %         plot3(x(k), y(k), z(k), 'ws--', 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r')
% % 
% % end