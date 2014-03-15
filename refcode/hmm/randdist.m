function v=randdist(row,col,val,prob)
len=row*col;

prob_cum=cumsum(prob(:));
prob_cum=prob_cum/max(prob_cum);
prob_cum=[0;
          prob_cum(1:end-1)];
m=rand(1,len);
for c=1:len
    m(c)=sum(m(c)>prob_cum);
end
v=reshape(val(m),row,col);
