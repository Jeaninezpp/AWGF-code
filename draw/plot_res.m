clc
clear all
% res = load('purity.txt');
% res = load('100leaves_nmi.txt');
% res = load('buaa_nmi.txt');
% res = load('caltech7_nmi.txt');
 res = load('mfeat_nmi.txt');
% res = res';
% res = res / 100;
grid on
plot((0.1:0.1:0.7),res(1,:),'b-','LineWidth', 1.5,'Marker','o','MarkerFaceColor','b');
hold on
plot((0.1:0.1:0.7),res(2,:),'m-','LineWidth', 1.5,'Marker','v','MarkerFaceColor','b');
plot((0.1:0.1:0.7),res(3,:),'g-','LineWidth', 1.5,'Marker','+','MarkerFaceColor','g');
plot((0.1:0.1:0.7),res(4,:),'k-','LineWidth', 1.5,'Marker','<','MarkerFaceColor','k');
plot((0.1:0.1:0.7),res(5,:),'r-','LineWidth', 1.5,'Marker','>','MarkerFaceColor','m');
% plot((5:5:30),res(6,:),'m--','LineWidth', 1.5,'Marker','h','MarkerFaceColor','m');
% 
% plot((5:5:30),res(7,:),'c--','LineWidth', 1.5,'Marker','p','MarkerFaceColor','c');
% plot((5:5:30),res(8,:),'c-','LineWidth', 1.5,'Marker','s','MarkerFaceColor','c');
% plot((5:5:30),res(9,:),'r-','LineWidth', 1.5,'Marker','o','MarkerFaceColor','r');

gca1 = legend({'DAIMC','IMG','MIC','INMF\_AGL','Ours'},'Location','southwest');
%set(gca,'YLim',[0.15 0.5]);
title('mfeat','FontSize',15);
ylabel('NMI','FontSize',15);
xlabel('Incomplete Ratio','FontSize',15);
set(gca,'PlotBoxAspectRatio',[1 1 1]);