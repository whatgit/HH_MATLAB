
cloud = load('X:\MT\Mapfilterstuff\output\pointcloud2.txt'); 

xx = cloud(:,1);
yy = cloud(:,2);
zz = cloud(:,3);
n=1;

for kk = 1:length(xx)
    if((abs(xx(kk))<100)&&((abs(yy(kk))<100))&&((zz(kk)<50)&&(zz(kk)>-50)))
%         if ((x(kk)<-0.4)||(abs(y(kk))>23)||(x(kk))>1)
        x(n,1) = xx(kk);
        y(n,1) = yy(kk);
        z(n,1) = zz(kk);
        data(n,:) = [x(n,1) y(n,1) z(n,1)];
        n = n+1;
%         end
    end

end


nn = 1;

for kk = 1:length(xx)
    if((abs(xx(kk))<5)&&((abs(yy(kk))<100))&&((zz(kk)<5)&&(zz(kk)>-2)))
%         if ((x(kk)<-0.4)||(abs(y(kk))>23)||(x(kk))>1)
        x1(n,1) = xx(kk);
        y1(n,1) = yy(kk);
        z1(n,1) = zz(kk);
        data1(n,:) = [x1(n,1) y1(n,1) z1(n,1)];
        n = n+1;
%         end
    end

end

% 
% 
% 
% figure;
% scatter(data(:,1),data(:,2))
% title('(a) Raw 3D Point Cloud')
% xlabel('X (m)'); 
% ylabel('Y (m)');
% zlabel('Z (m)');
% axis([-25 25 -10 40]);
% 
% figure;
% scatter(data(:,2),data(:,3))
% title('(b) Raw 3D Point Cloud')
% xlabel('Y (m)');
% ylabel('Z (m)');
% axis([-10 40 -40 10]);
% figure;
% % subplot(2,2,3);
% scatter(data1(:,1),data1(:,2))
% title('(c) Point Cloud without Outliers')
% xlabel('X (m)');
% ylabel('Y (m)');
% axis([-25 25 -10 40]);
% figure;
% subplot(2,2,4); scatter(data1(:,2),data1(:,3))
% title('(d) Point Cloud without Outliers')
% xlabel('Y (m)');
% ylabel('Z (m)');
% axis([-10 40 -40 10]); 

figure;
scatter(data1(:,1),data1(:,2))
xlabel('X (m)');
ylabel('Y (m)');
zlabel('Z (m)');
axis ([-5 5 -3 35]);

% figure;
% scatter(data1(:,2),data1(:,3))
% xlabel('Y (m)');
% ylabel('Z (m)');
% axis ([-3 35 -5 10]);


