function [Zstar,obj] = algorithm(X, G, V, invVV, Z, Zstar, Y, k, lmd1, lmd2, lmd3, max_iter)
% N: totall numbers
% k : clusters
num_views = length(X);
alpha = (1/num_views)*ones(1,num_views);

U= cell(1,num_views);

tol=1e-5;%convergence

for iter = 1:max_iter
    %disp(['Iter-', num2str(iter)]);
    sumofobj = 0;
    % ---------- Update Ui ----------%
    for iv = 1:num_views
        %U{iv}=X{iv}*V{iv}'*invVV{iv};
        Xtmp   = X{iv};
        Vtmp   = V{iv};
        invtmp = invVV{iv};
        
        Xtmp    = gpuArray(Xtmp);
        Vtmp    = gpuArray(Vtmp);
        invtmp  = gpuArray(invtmp);
        Utmp= Xtmp * Vtmp'*invtmp;
        Utmp   = gather(Utmp);
        
        U{iv}  = Utmp;
    end
    
    % ---------- Update Vi ----------%
    V=updateV(X,U,V,Z,lmd1,num_views);
    for iv = 1:num_views
        vtemp = V{iv};
        vtemp = gpuArray(vtemp);
        invtemp = inv(vtemp*vtemp');
        invtemp =gather(invtemp);
        invVV{iv} = invtemp;
    end
    % ---------- Update Zi ----------%
%>> QP solution
%     for iv =1:num_views
%         Znum = size(Z{iv},1);
%         
%         ZH = G{iv}*Zstar*G{iv}';
%         VV=V{iv}'*V{iv};
%         H = lmd1*VV+lmd2*alpha(iv)^2*eye(Znum);
%         f = lmd1*VV+lmd2*alpha(iv)*ZH;
%         
%         A=-1*eye(Znum);
%         Aeq = ones(1,Znum);
%         
%         Ztmp = Z{iv};
%         % Ztmp = gpuArray(Ztmp);  quadprog not support
% 
%         for col = 1:Znum
%             %X0=1/Znum*ones(Znum,1);
%             Ztmp(:,col)=quadprog(H,f(:,col)',A,zeros(Znum,1),Aeq,1);
%         end
%         % Ztmp = gather(Ztmp); 
%         Z{iv}=Ztmp;
%     end
%>>closed solution

    % ---------- Update Zi ----------%
    for iv = 1:num_views
        Znum = size(Z{iv},1);
        Zi = Z{iv};
        GZG = G{iv}*Zstar*G{iv}';

        Vi=V{iv};
        Vi = gpuArray(Vi);
        VV = Vi'*Vi;
        VV = gather(VV);

        Q = lmd1 * VV + lmd2 * (alpha(iv)^2) * eye(Znum);

        Q = gpuArray(Q);
        QQinv = inv(Q+Q');
        QQinv = gather(QQinv);

        for j = 1:Znum
            PP = 2 *  lmd1 * VV(:,j) + 2 * lmd2 * alpha(iv) * GZG(:,j);
            Zi(:,j) = QQinv * PP;
        end
        % constraint 
        Zi(Zi<0)=0;
        colsum = sum(Zi,1);
        colsum_diag = diag(colsum);
        Zi = Zi * colsum_diag^-1;
        % symmetry
        %Z{iv} = (Zi+Zi')/2;
        Z{iv} = Zi;
    end

    % ---------- Update alpha_i ----------%
    sumTrZHZZ = 0;
    suma = 0;
    multipTrZZ = 1;
    for iv = 1:num_views
        HH = G{iv}*Zstar*G{iv}';
        TrZH(iv) = trace(Z{iv}'*HH);
        TrZZ(iv) = trace(Z{iv}'*Z{iv});
        TrZHZZ = TrZH(iv)/TrZZ(iv);
        sumTrZHZZ = sumTrZHZZ +TrZHZZ;
        
        suma= suma+1/TrZZ(iv);
    end

    beta = (sumTrZHZZ-1)*2/suma;
    
    for iv =1:num_views
        alpha(iv) = (2*TrZH(iv)-beta)/(2*TrZZ(iv));
    end
    
    % ---------- Update Z* with for loooop ----------%
    % for each Z*_{kq}
%     for kk = 1:N
%         for qq = 1:N
%             for iv = 1:num_views
%             % find the kk qq sample in each view
%                 jj   = find(G{iv}(:,kk)==1);
%                 pp = find(G{iv}(:,qq)==1);
%                 if jj & pp
%                     rZ(iv) = Z{iv}(jj,pp);
%                 else
%                     rZ(iv) = 0;
%                 end
%             end
%             r = nnz(rZ);% number of non-zeros.
%             
%             sum_arZjp = 0;
%             if r == 0
%                 sum_arZjp = 0;
%             else
%                 for iv =1:r
%                     sum_arZjp = sum_arZjp + alpha(iv)*r*rZ(iv);
%                 end
%             end
% 
%             if sum_arZjp > lmd3/(2*lmd2)
%                 Zstar(kk,qq) = (2*lmd2*sum_arZjp - lmd3)/(2*lmd2*(r^2)+1e-10);
%             elseif sum_arZjp < -lmd3/(2*lmd2)
%                 Zstar(kk,qq) = (2*lmd2*sum_arZjp + lmd3)/(2*lmd2*(r^2)+1e-10);
%             else
%                 Zstar(kk,qq) = 0;
%             end
%             %disp(Zstar(kk,qq))
%         end
%     end
    
    % ---------- Update Z* without loooop ----------%
    sum_alpha_Zi = 0;
    sum_r = 0;
    for iv = 1:num_views
        bigZi{iv} = G{iv}'*Z{iv}*G{iv};
        sum_alpha_Zi = sum_alpha_Zi + alpha(iv) .* bigZi{iv};
        
        num_r{iv} = (bigZi{iv} ~= 0);
        sum_r = sum_r + num_r{iv};
    end
    Zstar( sum_r==0 ) = 0;
    A = (2 * lmd2 * sum_r .* sum_alpha_Zi - lmd3)./(2 * lmd2 * (sum_r.^2)+1e-10);
    B = (2 * lmd2 * sum_r .* sum_alpha_Zi + lmd3)./(2 * lmd2 * (sum_r.^2)+1e-10);
    Zstar( sum_alpha_Zi > lmd3 ./ (2 * lmd2 .* sum_r) ) = A( sum_alpha_Zi > lmd3 ./ (2 * lmd2 .* sum_r) );
    Zstar( sum_alpha_Zi < -lmd3 ./ (2 * lmd2 .* sum_r) ) = B( sum_alpha_Zi < -lmd3 ./ (2 * lmd2 .* sum_r) );

%     Zstar(Zstar<0)=0;
%     Zstar = (Zstar+Zstar')/2;
    
    Zstar  = (abs(Zstar)+abs(Zstar)')/2;

    %iter_clustering(iter,:) = spcclust(Zstar, k, Y);

    %====obj=====
    for iv = 1:num_views
        NMFterm = norm(X{iv}-U{iv}*V{iv},'fro');
        selfrepres = lmd1*norm(V{iv}-V{iv}*Z{iv},'fro');
        ZZloss      = lmd2*norm(alpha(iv)*Z{iv}-G{iv}*Zstar*G{iv}','fro');
        sumofview =  NMFterm + selfrepres + ZZloss;
        %fprintf('View-%d   NMFterm:%g \t selfrepres:%g \t ZZloss:%g \t  sumofview:%g \n',iv, NMFterm,selfrepres,ZZloss,sumofview);
        sumofobj = sumofobj + sumofview;
    end
    regularZ = lmd3*sum(sum(abs(Zstar)));
    obj(iter) = sumofobj + regularZ;
    %fprintf('regularZ:%g \t obj:%g \n',regularZ,obj(iter));

    if iter>19 && ((abs(obj(iter)-obj(iter-1))/obj(iter-1) < tol) || obj(iter)<=tol)
    %if iter == 20
        %fprintf('Objective value converge to %g at iteration %d before the maxIteration reached \n',obj(iter),iter);
        break;
    end
end
end
