file=input('Enter name of .txt file to take input from:\n','s');
np=fscanf(fopen(file),'%d',[1 1]);
data=readmatrix(file);
x=data(1:np,1);
y=data(1:np,2);
m=data(np+1,1);
x_=data(np+2:np+m+1,1);
choice=input('Choose method:\n1.Lagrange polynomials\n2.Natural Cubic Spline\n');

if(choice==1)
    l=cell(np,1);
    for i=1:np
        l{i}=@(z) 1;
        for k=1:np
            if(i~=k)
                l{i}=@(z) l{i}(z).*((z-x(k))/(x(i)-x(k)));
            end
        end
    end
    f=@(z) 0;
    for k=1:np
        f=@(z)f(z)+y(k).*l{k}(z);
    end

    lag=fopen('Lagrange.txt','w');
    fprintf(lag,'Lagrange Polynomials\n');
    for i=1:m
        fprintf(lag,'%f %f',x_(i),f(x_(i)));
        fprintf(lag,'\n');
    end
    fprintf('Output saved to file "Lagrange.txt"\n');
    hold on

    plot(x,y,'.','markersize',12);
    fplot(f,[x(1),x(np)]);
    
    hold off

end

if(choice==2)
    s_dp=zeros(np,1);
    s_dp(1)=0;
    s_dp(np)=0;
    vec=zeros(np-2,1);
    mat=zeros(np-2,np-2);
    mat(1,1)=2*(x(3)-x(1));
    mat(1,2)=(x(3)-x(2));
    mat(np-2,np-3)=(x(np-1)-x(np-2));
    mat(np-2,np-2)=2*(x(np-1)-x(np-3));

    for i=2:np-3
        mat(i,i-1)=x(i+1)-x(i);
        mat(i,i)=2*(x(i+2)-x(i));
        mat(i,i+1)=x(i+2)-x(i+1);
    end
    
    for i=1:np-2
        vec(i,1)=6*((y(i+2)-y(i+1))/(x(i+2)-x(i+1)))-6*((y(i+1)-y(i))/(x(i+1)-x(i)));
    end

    s_dp(2:np-1)=linsolve(mat,vec);
    
    spline=cell(np-1,1);
    for j=1:np-1
            a=@(z)(s_dp(j).*(x(j+1)-z).^3 + s_dp(j+1).*(z-x(j)).^3)/(6*(x(j+1)-x(j)));
            b=@(z)(y(j)/(x(j+1)-x(j)) - (x(j+1)-x(j))*s_dp(j)/6).*(x(j+1)-z);
            c=@(z)(y(j+1)/(x(j+1)-x(j)) - (x(j+1)-x(j))*s_dp(j+1)/6).*(z-x(j));
            spline{j}=@(z) a(z)+b(z)+c(z);
    end
    
    y_=zeros(m,1);
    for i=1:m
        for j=1:np-1
            if(x_(i)>x(j) && x_(i)<x(j+1))
                y_(i)=spline{j}(x_(i));
                break;
            end
        end
    end
    spl=fopen('Spline.txt','w');
    fprintf(spl,'Spline Interpolation\n');
    for i=1:m
        fprintf(spl,'%f %f',x_(i),y_(i));
        fprintf(spl,'\n');
    end
    fprintf('Output saved to file "Spline.txt"\n');

    hold on

    plot(x,y,'.','markersize',12);
    for i=1:np-1
        fplot(spline{i},[x(i),x(i+1)]);
    end

    hold off
end