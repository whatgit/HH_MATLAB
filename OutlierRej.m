function dataout = OutlierRej(dd)

% Function   : OutlierRej
% 
% Purpose    : Cut off enormous numbers and tags
% 
% Parameters : dd, Point cloud data
% 
% Return     : dataout, pointcloud data
%              

%    dd = removeoutliers(din);
   
   x = dd(:,1);
   y = dd(:,2);
   z = dd(:,3);
   n =1;
   
for kk = 1:length(x)
    if((abs(x(kk))<3)&&((abs(y(kk))<100))&&((z(kk)<5.5)&&(z(kk)>-1.2)))
        if ((x(kk)<-0.8)||(x(kk))>0.8) % Cut off the tags
        dataout(n,1) = x(kk);
        dataout(n,2) = y(kk);
        dataout(n,3) = z(kk);
        n = n+1;
        end
    end

end




end