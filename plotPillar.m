function NA = plotPillar(pillarL,pillarR,pillarcloud)
% Function   : plotPillar
% 
% Purpose    : Plot pillar in cylinder
% 
% Parameters : pillarL,pillarR: left and right pillar
% 
% Return     : None
%           

figure;
for k= 1:(length(pillarL))
[X,Y,Z]=cylinder1([pillarL(k,1) pillarL(k,2) 0],[pillarL(k,1) pillarL(k,2) 4],0.15);
aa = surf(X,Y,Z);
set(aa, 'edgecolor','none')
% hold on;
% set(aa,'FaceColor',[0 0 1],'FaceAlpha',0.8);
% camlight
% shading interp
% lighting gouraud
end

for k= 1:(length(pillarR))
[X,Y,Z]=cylinder1([pillarR(k,1) pillarR(k,2) 0],[pillarR(k,1) pillarR(k,2) 4],0.15);
aa = surf(X,Y,Z);
set(aa, 'edgecolor','none')
hold on;
% set(aa,'FaceColor',[0 0 1],'FaceAlpha',0.8);
% camlight
% shading interp
% lighting gouraud
end
% hold off;
% scatter3(pillarcloud(:,1),pillarcloud(:,2),pillarcloud(:,3),'.','g');hold on;

title('3D Pillar');
xlabel('X (m)');
ylabel('Y (m)');
axis ([-3 3 -1 32 0 3]);
% axis equal;

NA = 1;

end