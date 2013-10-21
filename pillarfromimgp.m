function [pl pr Corrwidth] = pillarfromimgp(leftp,rightp,rl,ll,pillarcloud)

    for kk = 1:length(leftp)
        pl(kk,:) = [ll(1,1) leftp(kk,2)];
    end

    for kk = 1:length(rightp)
        pr(kk,:) = [rl(1,1) rightp(kk,2)];
    end

%     NA = plotPillar(pl,pr,pillarcloud);

    Corrwidth = abs(mean(leftp(:,1)) - mean(rightp(:,1)));

figure;
scatter(pl(:,1),pl(:,2)); hold on;
scatter(pr(:,1),pr(:,2)); hold on;
title('2D Point Cloud');
xlabel('X (m)');
ylabel('Y (m)');
axis equal
end