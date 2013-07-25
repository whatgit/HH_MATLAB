%{
Function to plot histogram of data
%}
function histPlot()
   
data = load('D:\Thesis\MATLAB\result\Node Results\data_4tag.txt');
[row col] = size(data);
calculated_data = zeros(1,row);
for i=1:row
   calculated_data(1,i) = sqrt((data(i,2)*data(i,2))+(data(i,3)*data(i,3)));
   calculated_data(1,i)
end


% data_8tag_ccw = [%nav, nav_cmd, all
%                  [calculated_data(1,1), calculated_data(1,2), calculated_data(1,3)];  %1st tag
%                  [calculated_data(1,4), calculated_data(1,5), calculated_data(1,6)];  %2nd tag
%                  [calculated_data(1,7), calculated_data(1,8), calculated_data(1,9)];  %3rd tag
%                  [calculated_data(1,10), calculated_data(1,11), calculated_data(1,12)];  %4th tag
%                  [calculated_data(1,13), calculated_data(1,14), calculated_data(1,15)];  %5th tag
%                  [calculated_data(1,16), calculated_data(1,17), calculated_data(1,18)];  %6th tag
%                  [calculated_data(1,19), calculated_data(1,20), calculated_data(1,21)];  %7th tag
%                  [calculated_data(1,22), calculated_data(1,23), calculated_data(1,24)];  %8th tag
%                  ]
%   figure;
%   hold on;
% 
%  bar(data_8tag_ccw,0.75,'grouped');
%  set(gca, 'XTick',1:8, 'XTickLabel',{'Tag 1' 'Tag 2' 'Tag 3' 'Tag 4' 'Tag 5' 'Tag 6' 'Tag 7' 'Tag 8'})
%  legend('IMU','IMU & Control Command','IMU & Control Command & PTAM');
%  ylabel('Distance Error (meter)');
%  
 data_4tag_ccw = [%nav, nav_cmd, all
                 [calculated_data(1,1), calculated_data(1,2), calculated_data(1,3)];  %1st tag
                 [calculated_data(1,4), calculated_data(1,5), calculated_data(1,6)];  %2nd tag
                 [calculated_data(1,7), calculated_data(1,8), calculated_data(1,9)];  %3rd tag
                 [calculated_data(1,10), calculated_data(1,11), calculated_data(1,12)];  %4th tag
                 ]
  figure;
  hold on;

 bar(data_4tag_ccw,0.75,'grouped');
 set(gca, 'XTick',1:8, 'XTickLabel',{'Tag 1' 'Tag 2' 'Tag 3' 'Tag 4'})
 legend('IMU','IMU & Control Command','IMU & Control Command & PTAM');
 ylabel('Distance Error (meter)');

end

