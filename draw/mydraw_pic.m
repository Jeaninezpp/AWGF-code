clc;
clear;

data = struct2cell(load('allresult.mat'));

datanum = length(data);
% 100leaves 10-60 [1,2,3,4,5,6]
% buaa　      10-70 [7,8,9,10,11,12,13]
% caltech7   10-70 [14,15,16,17,18,19,20]
% mfeat       10-70 [21,22,23,24,25,26,27]
% orl            10-70 [28,29,30,31,32,33,34]
% orlRn        10-70 [35,36,37,38,39,40,41]


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
for ddi = 7:13
    if ddi>=1 & ddi <=13
        miss = ddi-6;
    elseif ddi>=14 & ddi <=20
        miss = ddi-13;
    elseif ddi>=21 & ddi <=27
        miss = ddi-20;
    end
        
    load(['./para/',num2str(ddi),'.mat']);
    
    showresults = acc(1:6,1:6)*100;
    x = [1000, 100, 10, 0.1, 0.01, 0.001];
    y = [0.001, 0.01, 0.1, 10, 100, 1000];
    bar3(showresults);
    xlabel('$\lambda1$','interpreter','latex','fontsize',13),ylabel('$\lambda2$','interpreter','latex','fontsize',15), zlabel('ACC(%)','fontsize',15);
    
    switch ddi
        case {1,2,3,4,5,6}
            titlename = '100Leaves';
        case {7,8,9,10,11,12,13}
            titlename = 'BUAA'; 
        case {14,15,16,17,18,19,20}
            titlename = 'Caltech7';
        case{21,22,23,24,25,26,27}
            titlename = 'mfeat';
    end
    
    all_title = [titlename, ' missing',num2str(miss * 10),'%'];
    title(all_title,'fontsize',13);
    set(gca,'xticklabel',x,'yticklabel',y);
    saveas(gcf, ['./para/newsenfig/',all_title,'_acc.eps'], 'psc2')
    close(gcf);
    
    %%%%%%%%%
    showresults = nmi(1:6,1:6)*100;
    x = [1000, 100, 10, 0.1, 0.01, 0.001];
    y = [0.001, 0.01, 0.1, 10, 100, 1000];
    bar3(showresults);
    xlabel('$\lambda1$','interpreter','latex','fontsize',13),ylabel('$\lambda2$','interpreter','latex','fontsize',15), zlabel('NMI(%)','fontsize',15);
    
    title(all_title,'fontsize',13);
    set(gca,'xticklabel',x,'yticklabel',y);
    saveas(gcf, ['./para/newsenfig/',all_title,'_nmi.eps'], 'psc2')
    close(gcf);
    clear acc nmi
end