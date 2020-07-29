clear
clc
warning off;

% path = 'E:\codes -access -hu\results';
% addpath(genpath(path));
dataname = 'WikipediaArticles';

% file = load([dataname,'.mat']);
% accval10 =file.accval10;
% % showresults = accval10(21:31,21:31)*100+10;
% showresults = accval10(5:15,5:15)*100+10;
% x = [-5:1:5];
% y = x;
% bar3(showresults);
% xlabel('$2^{\beta}$','interpreter','latex','fontsize',20),ylabel('$2^{\lambda}$','interpreter','latex','fontsize',20), zlabel('ACC(%)','fontsize',20);
% title('Flower17','fontsize',20);
% set(gca,'xticklabel',x,'yticklabel',y);


file = load([dataname,'.mat']);
acc =file.Nmi;
% showresults = accval10(21:31,21:31)*100+10;
%% 2-para
for i =1:6
    for j= 1:6
        acc(i,j)=A((6*i-6)+j,4);
        nmi(i,j)=A((6*i-6)+j,5);
    end
end
showresults = acc(1:6,1:6)*100;
x = [1000, 100, 10, 0.1, 0.01, 0.001];
y = [0.001, 0.01, 0.1, 10, 100, 1000];
bar3(showresults);
xlabel('$\lambda1$','interpreter','latex','fontsize',15),ylabel('$\lambda2$','interpreter','latex','fontsize',15), zlabel('ACC(%)','fontsize',15);
title('Caltech7','fontsize',15);
set(gca,'xticklabel',x,'yticklabel',y);

showresults = nmi(1:6,1:6)*100;
x = [1000, 100, 10, 0.1, 0.01, 0.001];
y = [0.001, 0.01, 0.1, 10, 100, 1000];
bar3(showresults);
xlabel('$\lambda2$','interpreter','latex','fontsize',15),ylabel('$\lambda3$','interpreter','latex','fontsize',15), zlabel('NMI(%)','fontsize',15);
title('One-complete_ 100leaves_ miss10%','fontsize',15);
set(gca,'xticklabel',x,'yticklabel',y);
%% 1-para
for k=1:6
    acc1para(k)=A(k,4);
end
showresults = acc1para(1:6)*100;
%x = [0.001, 0.01, 0.1, 10, 100, 1000];
x = [1000, 100, 10, 0.1, 0.01, 0.001];
bar(showresults);
xlabel('$\lambda2$','interpreter','latex','fontsize',12),ylabel('ACC(%)','fontsize',12);
title('One-complete_ 100leaves_ miss50%','fontsize',12);
set(gca,'xticklabel',x,);

showresults = nmi(1:6,1:6)*100;
x = [1000, 100, 10, 0.1, 0.01, 0.001];
y = [0.001, 0.01, 0.1, 10, 100, 1000];
bar3(showresults);
xlabel('$\lambda2$','interpreter','latex','fontsize',15),ylabel('$\lambda3$','interpreter','latex','fontsize',15), zlabel('NMI(%)','fontsize',15);
title('One-complete_ 100leaves_ miss10%','fontsize',15);
set(gca,'xticklabel',x,'yticklabel',y);
