% plot script
x = transpose( network(:,[2,4]));
y = transpose( network(:, [3, 5]));
hold on
plot(x,y,'Color',[1,0,0])
 scatter(joint(:,3),joint(:,4),20,joint(:,9),'filled')

 hold off