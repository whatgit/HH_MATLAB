function [leftp rightp] = estimatepillar(tagdata)

% Function   : estimatepillar
% 
% Purpose    : Calculate estimated pillar postion in global frame
% 
% Parameters : tagdata, pillar projection in image frame
% 
% Return     : leftp, pillar position on left side
%              rightp, pillar position on right side

% tagdata = load('X:\MT\Mapfilterstuff\output\warehouse4_pointcloud4.txt');

pos = [tagdata(:,2) tagdata(:,3)];
z = tagdata(:,4);
heading = tagdata(:,5);

%% Set Para : Coef for Cam -> Pixel
% Mcx = 0.871576; 
% Mcy = 0.493151;

Mcx = 16/sqrt(337);
Mcy = 9/sqrt(337);

%% Cal
ppxl1 = [tagdata(:,7) tagdata(:,8)];
ppxl2 = [tagdata(:,9) tagdata(:,10)];

% pwc1(:,1) =((tagdata(:,5)-320)/640).*(Mcx.*z.*2);
% pwc1(:,2) =(-(tagdata(:,6)-160)/320).*(Mcy.*z.*2);
% pwc2(:,1) =((tagdata(:,7)-320)/640).*(Mcx.*z.*2);
% pwc2(:,2) =(-(tagdata(:,8)-160)/320).*(Mcy.*z.*2);

pwc1(:,1) = ((tagdata(:,7)-320.0).*((2.0*tagdata(:,4))/734.3));
pwc1(:,2) = ((180.0-tagdata(:,8)).*((2.0*tagdata(:,4))/734.3));
pwc2(:,1) = ((tagdata(:,9)-320.0).*((2.0*tagdata(:,4))/734.3));
pwc2(:,2) = ((180-tagdata(:,10)).*((2.0*tagdata(:,4))/734.3));

% pwc1(:,1) =tagdata(:,2)+((tagdata(:,5)-320)/640).*(Mcx.*z.*2);
% pwc1(:,2) =tagdata(:,3)+((tagdata(:,6)-160)/320).*(Mcy.*z.*2);
% pwc2(:,1) =tagdata(:,2)+((tagdata(:,7)-320)/640).*(Mcx.*z.*2);
% pwc2(:,2) =tagdata(:,3)+((tagdata(:,8)-160)/320).*(Mcy.*z.*2);

% pos1 = [pwc1(:,1) pwc1(:,2)] * R;
% pos2 = [pwc2(:,1) pwc2(:,2)] * R;

%% Rotate 

iIndex = size(heading,1);

for i=1:iIndex
    alfa = heading(i,1)*3.14/180;
    R = [cos(alfa) -sin(alfa); sin(alfa) cos(alfa)];
%     pos1(i,:) = pwc1(i,:) * R;
%     pos2(i,:) = pwc2(i,:) * R;
    pos1(i,:) = pwc1(i,:) ;
    pos2(i,:) = pwc2(i,:) ;
    
end

%% Convert to World Frame : CoorWorld = DroneC + CamFPos
pwc1(:,1) =tagdata(:,2)+pos1(:,1);
pwc1(:,2) =tagdata(:,3)+pos1(:,2);
pwc2(:,1) =tagdata(:,2)+pos2(:,1);
pwc2(:,2) =tagdata(:,3)+pos2(:,2);
    
    
%% Plot Pillar 
    figure;
    
    lidx = 1;
    ridx = 1;
    for i = 1:iIndex

%     plot(pwc1(i,1),pwc1(i,2),'--rs','LineWidth',0.1);hold on;
%                    
%                    plot(pwc2(i,1),pwc2(i,2),'--rs','LineWidth',0.1); hold on;
%     
    plot([pwc1(i,1) pwc2(i,1)],[pwc1(i,2) pwc2(i,2)]);hold on;

    if( (pwc1(i,1) < pos(i,1)) && (pwc2(i,1) < pos(i,1)) )
       left(lidx,:) =  [pwc2(i,1) pwc2(i,2)];
       lidx = lidx +1;
    end
    
    if( (pwc1(i,1) > pos(i,1)) && (pwc2(i,1) > pos(i,1)) )
       right(ridx,:) =  [pwc1(i,1) pwc1(i,2)];
       ridx = ridx +1;
    end
        
    end
    
scatter(pos(:,1),pos(:,2),'r');
  scatter(left(:,1),left(:,2),'g'); hold on;
  scatter(right(:,1),right(:,2),'g');hold on;
axis equal
title('Pillar Projection From Image');
xlabel('X (m)');
ylabel('Y (m)');
    
%% Plot Trace
    
figure;
  scatter(left(:,1),left(:,2),'g'); hold on;
  scatter(right(:,1),right(:,2),'b');hold on;
 scatter(pos(:,1),pos(:,2),'r');
 title('Endpoint of Pillar Projection');
 xlabel('X (m)');
 ylabel('Y (m)');
 axis equal
 
idxss11 = 1;
for idxss = 1: length(right)
    if(abs(right(idxss,1))>0.3)
        withouttag(idxss11,:) = right(idxss,:);
        idxss11 = idxss11+1;
    end
    
end

right = [];
right = withouttag;
 
 %%
threshold = 1.2;
tmpx = 0;
tmpy = 0;
cnt = 1;
lidx = 1;
ridx = 1;
leftpillar = [];
rightpillar = [];


leftwidx = [];
nbr = 1;
leftwidx(1) = nbr

for i = 1:(length(left)-1)
 if( abs(left(i,2) - left(i+1,2)) >= threshold )
     nbr = nbr + 1;
     leftwidx(i+1) = nbr;
 end
 
 if( abs(left(i,2) - left(i+1,2)) < threshold )
     leftwidx(i+1) = nbr;
 end
 
end

rightwidx = [];
 nbr = 1;
rightwidx(1) = nbr
for i = 1:(length(right)-1)
 if( abs(right(i,2) - right(i+1,2)) >= threshold )
     nbr = nbr + 1;
     rightwidx(i+1) = nbr;
 end
 
 if( abs(right(i,2) - right(i+1,2)) < threshold )
     rightwidx(i+1) = nbr;
 end
 
end
 
n = 1; 

 for i = 1:(length(left))
     
     if leftwidx(i)== n
         tmpy = left(i,2) + tmpy;
         tmpx = left(i,1) + tmpx;
         cnt = cnt+1; 
     end
     
     if leftwidx(i) ~= n
        leftp(n,:) = [tmpx/cnt tmpy/cnt];
        cnt = 1;
        n = n+1;
         tmpy = left(i,2);
         tmpx = left(i,1);      
     end
    
     if i == (length(left))
         if leftwidx(i) == leftwidx(i-1)
             leftp(n,:) = [tmpx/cnt tmpy/cnt];         
         end  
     end
 end
 
 tmpx = 0;
tmpy = 0;
cnt = 1;
 
 
n = 1; 

 for i = 1:(length(right))
     
     if rightwidx(i)== n
         tmpy = right(i,2) + tmpy;
         tmpx = right(i,1) + tmpx;
         cnt = cnt+1; 
     end
     
     if rightwidx(i) ~= n
        rightp(n,:) = [tmpx/cnt tmpy/cnt];
        cnt = 1;
        n = n+1;
         tmpy = right(i,2);
         tmpx = right(i,1);      
     end
    
     if i == (length(right))
         if rightwidx(i) == rightwidx(i-1)
             rightp(n,:) = [tmpx/cnt tmpy/cnt];         
         end  
     end
 end
 
 
figure(1111);
scatter(leftp(:,1),leftp(:,2),'b','+'); hold on;
scatter(rightp(:,1),rightp(:,2),'b','+');hold on;  
scatter(pos(:,1),pos(:,2),'r','o'); hold on;
title('Pillar Estimated');
xlabel('X (m)');
ylabel('Y (m)');
axis equal


end