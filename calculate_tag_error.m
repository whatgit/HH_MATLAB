%{
Function to calculate the tag error
%}
function calculate_tag_error(type)

type = '8tags_ccw';
tag_scale = 15;

if strcmp(type,'4tags')
    tag_pos = zeros(4,2);
    tag_pos(1,1) = 0; %tag1_x
    tag_pos(1,2) = 0.7+0.1485; %tag1_y
    tag_pos(2,1) = 0; %tag2_x
    tag_pos(2,2) = 3.4+0.1485; %tag2_y
    tag_pos(3,1) = 0; %tag3_x
    tag_pos(3,2) = 6.1+0.1485; %tag3_y
    tag_pos(4,1) = 0; %tag4_x
    tag_pos(4,2) = 8.8+0.1485; %tag4_y
end

if strcmp(type,'8tags_cw')
    x_fixer = 4.62;
    y_fixer = 4.0;
    
    tag_pos = zeros(8,2);
    tag_pos(1,1) = mean([-1.025+x_fixer -1.025+(0.0210*tag_scale)+x_fixer]);    %tag1_x
    tag_pos(1,2) = mean([0.15+y_fixer 0.15+(0.0297*tag_scale)+y_fixer]);        %tag1_y
    tag_pos(2,1) = mean([0.935+x_fixer 0.935+(0.0210*tag_scale)+x_fixer]);      %tag2_x
    tag_pos(2,2) = mean([2.85+y_fixer 2.85+(0.0297*tag_scale)+y_fixer]);        %tag2_y
    tag_pos(3,1) = mean([0.265+x_fixer-(0.0297*tag_scale) 0.265+x_fixer]);      %tag3_x
    tag_pos(3,2) = mean([4.25+y_fixer 4.25+(0.0210*tag_scale)+y_fixer]);        %tag3_y
    tag_pos(4,1) = mean([-1.935+x_fixer-(0.0297*tag_scale) -1.935+x_fixer]);    %tag4_x
    tag_pos(4,2) = mean([4.74+y_fixer 4.74+(0.0210*tag_scale)+y_fixer]);        %tag4_y
    tag_pos(5,1) = mean([-4.505+x_fixer-(0.0210*tag_scale) -4.505+x_fixer]);    %tag5_x
    tag_pos(5,2) = mean([4.74+y_fixer-(0.0297*tag_scale) 4.74+y_fixer]);        %tag5_y
    tag_pos(6,1) = mean([-5.385+x_fixer-(0.0210*tag_scale) -5.385+x_fixer]);    %tag6_x
    tag_pos(6,2) = mean([3.28+y_fixer-(0.0297*tag_scale) 3.28+y_fixer]);        %tag6_y
    tag_pos(7,1) = mean([-3.835+x_fixer -3.835+(0.0210*tag_scale)+x_fixer]);    %tag7_x
    tag_pos(7,2) = mean([0.58+y_fixer-(0.0297*tag_scale) 0.58+y_fixer]);        %tag7_y
    tag_pos(8,1) = mean([-4.600+x_fixer-(0.0210*tag_scale) -4.600+x_fixer]);    %tag8_x
    tag_pos(8,2) = mean([-2.12+y_fixer-(0.0297*tag_scale) -2.12+y_fixer]);      %tag8_y
    
end

if strcmp(type,'8tags_ccw')
    tag_pos = zeros(8,2);
    tag_pos(1,1) = mean([-1.025 -1.025+(0.0210*tag_scale)]);    %tag1_x
    tag_pos(1,2) = mean([0.15 0.15+(0.0297*tag_scale)]);        %tag1_y
    tag_pos(2,1) = mean([0.935 0.935+(0.0210*tag_scale)]);      %tag2_x
    tag_pos(2,2) = mean([2.85 2.85+(0.0297*tag_scale)]);        %tag2_y
    tag_pos(3,1) = mean([0.265-(0.0297*tag_scale) 0.265]);      %tag3_x
    tag_pos(3,2) = mean([4.25 4.25+(0.0210*tag_scale)]);        %tag3_y
    tag_pos(4,1) = mean([-1.935-(0.0297*tag_scale) -1.935]);    %tag4_x
    tag_pos(4,2) = mean([4.74 4.74+(0.0210*tag_scale)]);        %tag4_y
    tag_pos(5,1) = mean([-4.505-(0.0210*tag_scale) -4.505]);    %tag5_x
    tag_pos(5,2) = mean([4.74-(0.0297*tag_scale) 4.74]);        %tag5_y
    tag_pos(6,1) = mean([-5.385-(0.0210*tag_scale) -5.385]);    %tag6_x
    tag_pos(6,2) = mean([3.28-(0.0297*tag_scale) 3.28]);        %tag6_y
    tag_pos(7,1) = mean([-3.835 -3.835+(0.0210*tag_scale)]);    %tag7_x
    tag_pos(7,2) = mean([0.58-(0.0297*tag_scale) 0.58]);        %tag7_y
    tag_pos(8,1) = mean([-4.600-(0.0210*tag_scale) -4.600]);    %tag8_x
    tag_pos(8,2) = mean([-2.12-(0.0297*tag_scale) -2.12]);      %tag8_y
end

data = load('D:\Thesis\MATLAB\matlab\data\edited_ccw_nav.txt');
[row col] = size(data);
for i=1:row
    disp(['Data No. : ' num2str(i) ...
        ', Tag No. : ' num2str(data(i,2)) ...
        ', Error in x : ' num2str(tag_pos(data(i,2),1)-data(i,3))...
        ', Error in y : ' num2str(tag_pos(data(i,2),2)-data(i,4))]);
end

end