function C = lspoly(x,y,M)

n = length(x);
F = zeros(n,M+1);
for k = 1:M+1
   F(:,k) = x'.^(k-1);
end
A = F'*F;
B = F'*y';
C = A\B;