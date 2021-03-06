% Pei Zhang 2020-07-29
% sum ||X-UV||+||V-VZ||+||aWZW-Z'||+b||Z'||_*
clear;
clc

resultdir = 'results/pairwise/one_parameter/';
if (~exist('results/pairwise/one_parameter','file'))
    mkdir('results/pairwise/one_parameter');
    addpath(genpath('results/pairwise/one_parameter/'));
end
%% Example Title
datadir = './setting/pairwised-data/';
dataname = {'buaaRnSp','caltech7','100Leaves','mfeatRnSp','ORL'}; 
datanum = length(dataname);
for datai = 1:datanum
    datafile = [datadir, cell2mat(dataname(datai))];
    load(datafile); % truth,data,per10-per70,
    fprintf('%s...\n',datafile);
    
    num_clusters = length(unique(truth));
    num_views    = length(data);
    num_sample  = length(truth);

    for per_in = 1:length(per)  % per incomplete ratio 
        in_ratio = (per_in*2-1)*10;
        percent = per{per_in};%percent = cell(1,10);
        disp(['missing',num2str(in_ratio),'%']);
        for folds = 1:1
            foldspath = [resultdir, char(dataname(datai)),'/Fold',num2str(folds)];
            if (~exist(foldspath,'file'))
                mkdir(foldspath);
                addpath(genpath([foldspath,'/']));
            end 
            
            index = percent{folds};
            savetxt = [resultdir ,'pairwise_',char(dataname(datai)),'_missing',num2str(in_ratio),'%','.txt'];
            dlmwrite(savetxt, ['Folds = ', num2str(folds)],'-append','delimiter','\t','newline','pc');
            %% get the incomplete Xi, Gi, Vi, Zi
            for iv = 1:num_views
                X_complete = data{iv};% di × n
                X_complete = NormalizeFea(X_complete,1);
                exist_index{iv}  = find(index(:,iv) == 1);
                X_incomplete{iv} = X_complete(:,exist_index{iv});% d_i * n_i
                ni_num(iv) = size(X_incomplete{iv},2);% num of existing sample

                % ===Generate n * n_i incomplete indicator: Gi===
                % If the i-th sample in n is the j-th sample in n_i, the G_ij =1.
                Gtmp = zeros(num_sample,length(exist_index{iv}));% n * n_i
                for gi = 1:length(exist_index{iv})
                    Gtmp(exist_index{iv}(gi),gi)=1;
                end
                G{iv}=Gtmp'; % n_i * n
            end
            clear X X_complete
            X=X_incomplete; % X \in di * ni
            clear X_incomplete

            Zcount = zeros(num_sample);
            Zbigsum = zeros(num_sample);
            for iv = 1:num_views
                % ===Initial Vi=== k*ni
                Vtmp = litekmeans(X{iv}',num_clusters,'MaxIter',100);% Vtmp: ni*1
                tmp = zeros(ni_num(iv),num_clusters);
                tmp(sub2ind(size(tmp),[1:ni_num(iv)],Vtmp'))=1;% ni*1->ni*k
                V{iv} = tmp';

                % ===Initial Zi===
                options = [];
                options.NeighborMode = 'KNN';
                options.k = 3;
                % options.WeightMode = 'Binary';
                Z1 = constructW(X{iv}',options);% input N*d dim X
                Z_ini{iv} = full(Z1);
                clear Z1

                % ===Initial Z*===
                Zbig{iv} = G{iv}'*Z_ini{iv}*G{iv};
                Zbig10    = Zbig{iv};
                Zbig10(Zbig10~=0) = 1;
                Zcount    = Zcount + Zbig10;
                Zbigsum  = Zbigsum + Zbig{iv};
            end
            Zstar_ini = Zbigsum./(Zcount+1e-10);
                % === finish initial Z*===
            clear  Gtmp Vtmp tmp


            % compute inverse GPU vers
            for iv = 1:num_views
                vtemp = V{iv};
                vtemp = gpuArray(vtemp);
                invtemp = inv(vtemp*vtemp');
                invtemp =gather(invtemp);
                invVV{iv} = invtemp;
            end

            % CPU vers
            % tic
            % for iv = 1:num_views
            %     vtemp = V{iv};
            %     vtemp = gpuArray(vtemp);
            %     invtemp = inv(vtemp*vtemp');
            %     invVV{iv}=gather(invtemp);
            % end
            % t2 = toc

            %% Method
            %
            %
            % V:  k*n_i
            % 
            % 
            max_iter =100;
            
            lambda1 = 1;
            lambda2=[1e3,1e2,1e1,1e-1,1e-2,1e-3];
            lambda3=[1e-3,1e-2,1e-1,1e1,1e2,1e3];
            lambda2=1;
            
            for i=1:length(lambda1)
                for j = 1:length(lambda2) 
                    for d = 1:length(lambda3)
                        for repi = 1:1
                            % disp(['Repeat:',num2str(repi)]);
                            tic;
                            [Zstar, obj]  = algorithm_v2(X, G, V, invVV, Z_ini, Zstar_ini, truth ,num_clusters, lambda1(i), lambda2(j), lambda3(d), max_iter);
                            metric(repi,:) = spcclust(Zstar, num_clusters, truth);
                            %fprintf('acc: %f \t nmi: %f \t Fscore: %f \t AR: %f \n', metric(repi,1),metric(repi,2),metric(repi,3),metric(repi,4));
                            one_repi_time(repi) = toc;
                            %disp(['one_repi_time:',num2str(one_repi_time(repi))]);
                        end
                        mean_one_repi_time = mean(one_repi_time);
                        Final_results=[];
                        Final_results(1) = mean(metric(:,1));% ACC
                        Final_results(2) = mean(metric(:,2));% NMI
                        Final_results(3) = mean(metric(:,3));% Fscore
                        Final_results(4) = mean(metric(:,4));% AR
                        %

                        disp(['lmd1: ',num2str(lambda1(i)),'    lmd2: ',num2str(lambda2(j)),'    lmd3: ',num2str(lambda3(d)),'    Acc: ',num2str(Final_results(1)),'    Nmi: ',num2str(Final_results(2)),'    One_repi_time:',num2str(one_repi_time(repi))]);
                        Final_results = [lambda1(i), lambda2(j), lambda3(d),Final_results];

                        % savetxt = [resultdir ,char(dataname),'.txt'];
                        dlmwrite(savetxt, Final_results ,'-append','delimiter','\t','newline','pc');
                        matname = [resultdir, char(dataname(datai)),'/Fold',num2str(folds),'/',num2str(lambda1(i)),'_',num2str(lambda2(j)),'_', num2str(lambda3(d)),'_',num2str(in_ratio),'%_.mat'];
                        save(matname, 'Final_results','Zstar','obj','mean_one_repi_time');
            %             savexls = [resultdir ,char(dataname),'.xls'];
            %             xlswrite(savexls, Final_results);
                    end
                end
            end
        end % folds end
    end % missing ratio end
end
