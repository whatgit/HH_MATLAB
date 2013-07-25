%{
plotResult function plot the trace and node position 


%}
function plotResult(type)

type = 'plotNodeTag';

figure;
xlabel('X (meter)');
ylabel('Y (meter)');
grid on, hold on;

if strcmp(type,'plotNodeTag')
    
    %Load the drone position data
    pose_data = load('D:\Thesis\MATLAB\matlab\data\set2\viso2_pose_dataset3.txt');
    [data_row data_col] = size(pose_data);

    plot(pose_data(:,5),-1*pose_data(:,6),'-xr','LineWidth',1);
    
    %Load the tag position data
    tag_data = load('D:\Thesis\MATLAB\matlab\data\set2\viso2_pose_dataset1.txt');
            plot(tag_data(:,3),tag_data(:,4),'--gs','LineWidth',2,...
                                'MarkerEdgeColor','k',...
                                'MarkerFaceColor','k',...
                                'MarkerSize',10);

    legend('PTAM & IMU & Control command', 'Tag Detected');

end

if strcmp(type,'plotTrace')
    
    %Plot it in green line only nav
    data = load('D:\Thesis\MATLAB\matlab\data\set2\dataset10\pose_nav.txt');
    [data_row data_col] = size(data);
    
    plot(data(:,4),data(:,5),'-xg','LineWidth',1);
    
    %Plot it in red line nav+cmd
    data = load('D:\Thesis\MATLAB\matlab\data\set2\dataset10\pose_nav_cmd.txt');
    [data_row data_col] = size(data);
    
    plot(data(:,4),data(:,5),'-xr','LineWidth',1);
    
    %Plot it in black line PTAM+nav+cmd
    data = load('D:\Thesis\MATLAB\matlab\data\set2\dataset10\pose_all.txt');
    [data_row data_col] = size(data);
    plot(data(:,4),data(:,5),'-xk','LineWidth',1);
    
    if withNewFigure == 1
        legend('IMU','IMU & Control command','PTAM & IMU & Control command');
    end

end

%create map type to plot on
mapCreate('start_tag_21');

%Configure axis if needed . . . some suggestion below
%axis([-4 4 -1 13]);    %axis for 4tags
%axis([-8 4 -4 8]);     %axis for 8tag_ccw
%axis([-7 20 -5 42]);  %axis for start_tag_1
%axis([-7 2 -5 42]);     %axis for start_tag_22
%axis([-7 2 -5 42]);   %axis for start_tag_21
%axis([-10 18 -5 42]);   %axis for start_tag_28
%axis([-15 2 -5 42]);  %axis for start_tag_16

end