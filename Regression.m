file=input('Enter name of .txt file to take input from:\n','s');
np=fscanf(fopen(file),'%d',[1 1]);
all=readmatrix(file);
x=all(:,1);
y=all(:,2);

m=input('Enter degree of polynomial:\n');
A=zeros(m+1,m+1);
B=zeros(m+1,1);

for i=1:m+1
    for j=1:m+1
        for k=1:np
            A(i,j)=A(i,j)+((x(k))^(i+j-2));
        end
    end

    for k=1:np
        B(i)=B(i)+y(k)*((x(k))^(i-1));
    end
end

C=linsolve(A,B);
f=@(x)0;
for i=1:m+1
    f=@(x)f(x)+C(i)*(x.^(i-1));
end

sr=0;
st=0;
mean=sum(y)/np;

for i=1:np
    sr=sr+(y(i)-f(x(i)))^2;
    st=st+(y(i)-mean)^2;
end

r_sq=(st-sr)/st;

reg=fopen('Regression.txt','w');
fprintf(reg,'Degree of polynomial = %d\nCoffeicients: ',m);
for i=1:m+1
    fprintf(reg,'%f ',C(i));
end
fprintf(reg,'\nR-sq = %f',r_sq);
fprintf('Output saved to file "Regression.txt"\n');

hold on

plot(x,y,'.','markersize',12);
fplot(f,[min(x(:)),max(x(:))]);

hold off




