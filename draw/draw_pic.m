clc;
clear;

data = struct2cell(load('../allresult.mat'));

datanum = length(data);
% 100leaves 10-60
% buaa　      10-70
% caltech7   10-70
% mfeat       10-70
% orl            10-70
% orlRn        10-70


%% convert txt into mat
% for datai = 1:datanum
%     now = data{datai};
%     for i =1:6
%         for j= 1:6
%             acc(i,j)=now((6*i-6)+j,4);
%             nmi(i,j)=now((6*i-6)+j,5);
%             save(['./para/',num2str(datai),'.mat'],'acc','nmi');
%         end
%     end
% end

%% sensitivity
for ddi = 33:41
    load(['./para/',num2str(ddi),'.mat']);
    
    showresults = acc(1:6,1:6)*100;
    x = [1000, 100, 10, 0.1, 0.01, 0.001];
    y = [0.001, 0.01, 0.1, 10, 100, 1000];
    bar3(showresults);
    xlabel('$\lambda2$','interpreter','latex','fontsize',13),ylabel('$\lambda3$','interpreter','latex','fontsize',15), zlabel('ACC(%)','fontsize',15);
    %title('One-complete_ 100leaves_ miss10%','fontsize',13);
    set(gca,'xticklabel',x,'yticklabel',y);
    saveas(gcf, ['./para/senfig/',num2str(ddi),'_acc.eps'], 'psc2')
    close(gcf);
    
    
    showresults = nmi(1:6,1:6)*100;
    x = [1000, 100, 10, 0.1, 0.01, 0.001];
    y = [0.001, 0.01, 0.1, 10, 100, 1000];
    bar3(showresults);
    xlabel('$\lambda2$','interpreter','latex','fontsize',13),ylabel('$\lambda3$','interpreter','latex','fontsize',15), zlabel('NMI(%)','fontsize',15);
    %title('One-complete_ 100leaves_ miss10%','fontsize',13);
    set(gca,'xticklabel',x,'yticklabel',y);
    saveas(gcf, ['./para/senfig/',num2str(ddi),'_nmi.eps'], 'psc2')
    close(gcf);
    clear acc nmi
end