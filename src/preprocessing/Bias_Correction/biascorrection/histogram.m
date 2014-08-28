function[h]=histogram(data)

data=data(:);
tam=length(data);
m=max(data)+1;
h=zeros(1,m);

for i=1:tam,
    f=floor(data(i));    
    if(f>0 & data(i)<=m)        
        a2=data(i)-f;
        a1=1-a2;
        h(f)  =h(f)  + a1;      
        h(f+1)=h(f+1)+ a2;                          
    end;
end;

h=conv(h,[1,2,3,2,1]);
h=h(3:(length(h)-2));

h=h/sum(h);


