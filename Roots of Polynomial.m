degree=input('Degree of the polynomial= ');
p=zeros(degree+1,1);
for i=1:degree+1
    fprintf('Coefficient of x^%d',i-1);
    p(i)=input('= ');
end
choice = input('Choose the number corresponding to method:\n1.Muller Method\n2.Bairstow Method\n');
max_err = input('Enter the maximum permissible error%\n');
max_it = input('Enter the maximum number of iterations to be performed\n');

if(choice==1)
    x=zeros(max_it+3);
    x(1)=input("Enter smallest value: ");
    x(2)=input("Enter middle value: ");
    x(3)=input("Enter largest value: ");
    err=zeros(1,max_it);
    iter=zeros(1,max_it);
    p=flip(p,1);
    f=@(a)polyval(p,a);
    for i=3:max_it+3
        d1=f(x(i))-f(x(i-1));
        d2=f(x(i-1))-f(x(i-2));

        num = (d1/(x(i)-x(i-1))) - (d2/(x(i-1)-x(i-2)));
        a=num/(x(i)-x(i-2));

        b = (f(x(i))-f(x(i-1)))/(x(i)-x(i-1)) + a*(x(i)-x(i-1));
        c = f(x(i));
        det=b^2-4*a*c;
        del_x1= -2*c/(b+sqrt(det));
        del_x2= -2*c/(b-sqrt(det));
        if(abs(del_x1)<abs(del_x2))
            x(i+1)=x(i)+del_x1;
        else
            x(i+1)=x(i)+del_x2;
        end
        
        if(i==3)
            continue
        end
        err(i)=abs((x(i+1)-x(i))/x(i+1))*100;
        iter(i)=i-3;
        if(err(i)<max_err)
            break
        end
    end
    fplot(f);
    fprintf('The root is %.10f +- %.10f %c\n',x(i+1),err(i),37);
end

if(choice==2)
    fprintf('Considering quadratic factor to be x^2-a1x-a0\n');
    a1=input('a1= ');
    a0=input('a0= ');
    
    n=degree+1;
    d=zeros(n,1);
    del=zeros(n-1,1);
    err=zeros(n,2);
    f=@(a)polyval(p,a);

    for i=1:max_it
        d(n)=p(n);
        d(n-1) = p(n-1) + a1*d(n);
        d(n-2)=p(n-2)+a1*d(n-1) + a0*d(n);
        del(n-1)=d(n);
        del(n-2)=d(n-1)+a1*del(n-1);
        for j=n-3:-1:1
           d(j)=p(j)+a1*d(j+1)+a0*d(j+2);
           del(j)=d(j+1)+a1*del(j+1)+a0*del(j+2);
           
        end
        
        A=[del(2),del(1); del(3),del(2)];
        B=[-d(1);-d(2)];
        delta = linsolve(A,B);
        a0=a0+delta(1);
        a1=a1+delta(2);

        err(i,1)=abs(d(1))*100;
        err(i,2)=abs(d(2))*100;

        if(max(err(i,1),err(i,2))<max_err) 
            break;
        end
        disp(d);
        disp(del);
    end
    
    det_sq=(a1*a1)+4*a0;
    r_1=0.5*(a1+sqrt(det_sq));
    r_2=0.5*(a1-sqrt(det_sq));

    if(det_sq<0)
        if(imag(r_1)<0)
                fprintf('The roots are %f%fi and %f+%fi\n',real(r_1),imag(r_1),real(r_2),imag(r_2));
        else
                fprintf('The roots are %f+%fi and %f%fi\n',real(r_1),imag(r_1),real(r_2),imag(r_2));
        end
    else
        fprintf('The roots are %f and %f\n',r_1,r_2);
    end
    fplot(f);

end
