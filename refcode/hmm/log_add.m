function Lc=log_add(La,Lb)

if (Lb>La)
    Lc=log_add(Lb,La);
    return;
end

if La==Lb
    Lc=La+log(2);
    return;
end

% iteration
tol=1.0e-10;

Lba=Lb-La;
f=exp(Lba);
if f<0.1
    Lc=La+f;
    n=2;
    s=-1;
    while 1
        tmp=exp(n*Lba)/n;
        Lc=Lc+s*tmp;
        if tmp<tol
            break;
        end
        n=n+1;
        s=-s;
    end
    %fprintf(1, 'method 1, n=%d\n',n-2);
    return    
else
    g=f/(2+f);
    Lc1=g;
    n=3;
    while 1
        tmp=g^n/n;
        Lc1=Lc1+tmp;
        if tmp<tol
            break;
        end
        n=n+2;
    end
    Lc=La+Lc1+Lc1;
    %fprintf(1, 'method 2, n=%d\n',(n-3)/2);
    return;
end

