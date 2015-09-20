clear all;
N=10;

k=zeros(N,1);
x=k;

for i=1:N
    k(i)=i-1;
    x(i)=-cos(pi/N.*(k(i)+0.5));
end

for i=1:N
    M(:,i)=x.^(i-1);
end


c0=zeros(N,1);
c1=c0;
c0(1)=1;
c1(2)=1;

C(:,1)=c0;
C(:,2)=c1;

for j=3:N
    for k=0:N-1
        C(k+1,j)=2*C(mod(k-1,N)+1,j-1)-C(k+1,j-2);
    end
end

C=sqrt(2)*C;
C(:,1)=C(:,1)/sqrt(2);
E=M*C


for i=1:N
    for j=1:N
        DOT(i,j)=dot(E(:,i),E(:,j))/N;
    end
end

x
M
C
DOT