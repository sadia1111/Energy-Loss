function y = horner(a,c)

n=length(a)-1;
y=a(n+1);

for k = n:-1:1
    y = a(k)+ c*y;
end