function Ls=log_sum(Lx)

[row col]=size(Lx);

% case of scalar
if (row==1) && (col==1)
    Ls=Lx;
    return;
end

% case of row vector
if row==1
    Ls=log_sum(Lx(:));
    return
end

% case of column vector
if col==1
    if row==2
        Ls=log_add(Lx(1),Lx(2));
        return;
    else
        Lx1=Lx(1:floor(row/2));
        Lx2=Lx(floor(row/2)+1:end);
        Ls=log_add(log_sum(Lx1),log_sum(Lx2));
        return;
    end
end

% case of matrix
Ls=zeros(1,col);
for c=1:col
    Ls(c)=log_sum(Lx(:,c));
end

return;
