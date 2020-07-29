function V=updateV(X,U,V,Z,lmd1,nv)
%
times = 0;
lastobj=0;
while 1
    times = times +1;
    sumVPminus=0;
    sumVPplus=0;
    sumUUVminus=0;
    sumUUVplus=0;
    sumUXminus=0;
    sumUXplus=0;
    for iv = 1:nv
        Ztmp = Z{iv};
        P = eye(size(X{iv},2))-Ztmp-Ztmp'+Ztmp*Ztmp';
        Pplus = (abs(P) + P)/2;
        Pminus = (abs(P) - P)/2;
        
        VPminus = lmd1*V{iv}*Pminus;
        VPplus    = lmd1*V{iv}*Pplus;
        
        UU = U{iv}'*U{iv};
        UUplus = (abs(UU)+UU)/2;
        UUminus = (abs(UU)-UU)/2;
        UUVplus    = UUplus * V{iv};
        UUVminus = UUminus * V{iv};
        
        UX = U{iv}'*X{iv};
        UXplus = (abs(UX)+UX)/2;
        UXminus = (abs(UX)-UX)/2;
        
        V{iv}=V{iv}.*sqrt((VPminus + UUVminus + UXplus)./(max(VPplus + UUVplus + UXminus, 1e-10)));
    end     
    

    
    obj = 0;
    for iv = 1:nv
        tmp = X{iv}-U{iv}*V{iv};
        obj = obj + sum(sum(tmp.^2));
    end
    if abs((obj-lastobj)/lastobj)<1e-4 || abs(obj-lastobj)>1e100 || times == 30
        break;
    end
    lastobj = obj;
end
end

