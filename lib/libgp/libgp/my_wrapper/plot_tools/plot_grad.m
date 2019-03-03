clear; clc

d = load('../data/grad1.dat');

close all

plot(d(:,1),d(:,2),'LineWidth',3);
hold on
plot(d(:,1),d(:,3),'--','LineWidth',3);

l=legend('  gp grad','  finite difference');
l.Location = 'best';


set(gca,'FontSize',20)

grid on
axis tight