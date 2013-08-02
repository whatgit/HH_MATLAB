%{
plot function for gui
%}
function gui_plot(pose_file,tag_file,plot_type)
    
    %Load the drone position data
%     pose_data = load('D:\Thesis\MATLAB\matlab\data\set2\outside3_pose_world.txt');
%     pose_data = load('D:\Thesis\MATLAB\matlab\data\set2\outside3_pose_world.txt');
    if strcmp(plot_type,'radio_plot_vo') || strcmp(plot_type,'radio_plot_both')
    pose_data = load(pose_file);
    [data_row data_col] = size(pose_data);
%     plot(pose_data(1,5),pose_data(1,6),'-gs','LineWidth',10);
%     plot(pose_data(data_row,5),pose_data(data_row,6),'-gs','LineWidth',10);
    plot(pose_data(:,4),pose_data(:,5),'-xr','LineWidth',1);
    
%     diff_x = pose_data(1,5)-pose_data(data_row,5)
%     diff_y = pose_data(1,6)-pose_data(data_row,6)
    end
  
    if strcmp(plot_type,'radio_plot_tag') || strcmp(plot_type,'radio_plot_both')
    %Load the tag position data
    tag_data = load(tag_file);
            plot(tag_data(:,3),tag_data(:,4),'--gs','LineWidth',2,...
                                'MarkerEdgeColor','k',...
                                'MarkerFaceColor','k',...
                                'MarkerSize',10);
    end
%     legend('PTAM & IMU & Control command', 'Tag Detected');

%create map type to plot on
% mapCreate('start_tag_21');

%Configure axis if needed . . . some suggestion below
%axis([-4 4 -1 13]);    %axis for 4tags
%axis([-8 4 -4 8]);     %axis for 8tag_ccw
%axis([-7 20 -5 42]);  %axis for start_tag_1
%axis([-7 2 -5 42]);     %axis for start_tag_22
%axis([-7 2 -5 42]);   %axis for start_tag_21
%axis([-10 18 -5 42]);   %axis for start_tag_28
%axis([-15 2 -5 42]);  %axis for start_tag_16

end