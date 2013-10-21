function [MeanDleft MeanDright varl varr]=calmeanvar(meanl,meanr)

% Function   : calmeanvar
% 
% Purpose    : Calculate mean value and variance
% 
% Parameters : meanl,meanr: left and right pillar
% 
% Return     : MeanDleft and MeanDright, Means
%              varl and varr, variance
%           


meanl;
meanr;
leftd = [];
rightd = [];

for l = 1:length(meanl)-1
    leftd(l,:) = [meanl(l+1,1)-meanl(l,1) meanl(l+1,2)-meanl(l,2)];
    
end

MaxDleft = [max(leftd(:,1)) max(leftd(:,2))]
MeanDleft = [mean(leftd(:,1)) mean(leftd(:,2))]

varl = var(leftd)
for r = 1:length(meanr)-1
    rightd(r,:) = [meanr(r+1,1)-meanr(r,1) meanr(r+1,2)-meanr(r,2)];
    
end

MaxDright = [max(rightd(:,1)) max(rightd(:,2))]
MeanDright= [mean(rightd(:,1)) mean(rightd(:,2))]
varr = var(rightd)

end
