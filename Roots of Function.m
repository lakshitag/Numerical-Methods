str = input('Enter your function f(x): ','s');
str = append('@(x) ',str);
f = str2func(str);
choice = input('Choose method:\na.Bisection\nb.False-position\nc.Fixed-Point\nd.Newton-Raphson\ne.Secant\n','s');
max_err = input('Enter the maximum permissible error%\n');
max_it = input('Enter the maximum number of iterations to be performed\n');

if(choice=='a')
    a = input("Enter smaller value: ");
    b = input("Enter larger value: ");  
    m=(a+b)/2;
    err=zeros(1,max_it);
    iter=zeros(1,max_it);
    for i=1:max_it
        iter(i)=i;
        if(f(m)==0) 
            break;
        end
        if(f(a)*f(m)<0) 
            b=m;
        end
        if(f(b)*f(m)<0)
            a=m;
        end
        m_prev=m;
        m=(a+b)/2;
        err(i)=(abs(m-m_prev))*100/abs(m);
        if (err(i)<max_err)
            break;
           
        end
        

    end
    line(iter(1:i),err(1:i),'LineWidth',1.75);
    fplot(f);
    fprintf('The root of %s = 0 is %.10f +- %.10f %c\n',extractAfter(func2str(f),')'),m,err(i),37);
end

if(choice=='b')
    a = input("Enter smaller value: ");
    b = input("Enter larger value: ");
    m = (a*f(b) - b*f(a))/(f(b)-f(a));
    err=zeros(1,max_it);
    iter=zeros(1,max_it);
    for i=1:max_it
        
        if(f(a)*f(m)<0)
            b=m;
        end
        if(f(b)*f(m)<0)
            a=m;
        end
        if(f(m)==0)
            break;
        end
        m_prev=m;
        m = (a*f(b) - b*f(a))/(f(b)-f(a));
        err(i)=abs(m-m_prev)*100/abs(m);
        iter(i)=i;
        if(err(i)<max_err)
            break;
        end
    end
    line(iter(1:i),err(1:i),'LineWidth',1.75);
    fplot(f);
    fprintf('The root of %s = 0 is %.10f +- %.10f %c\n',extractAfter(func2str(f),')'),m,err(i),37);
end
if(choice=='c')
    str_g=input('Enter g(x) such that f(x)=0 implies g(x)=x: ','s');
    str_g=append('@(x)',str_g);
    g=str2func(str_g);
    a = input('Enter initial guess: ');
    err=zeros(1,max_it);
    iter=zeros(1,max_it);
    for i=1:max_it
        err(i)=abs((g(a)-a)*100)/abs(g(a));
        iter(i)=i;
        if(err(i)<max_err)
            break;
        end
        a=g(a);
    end
    line(iter(1:i),err(1:i),'LineWidth',1.75);
    fplot(f);
    fprintf('The root of %s = 0 , that is, %s = x is %.10f +- %.10f %c\n',extractAfter(func2str(f),')'),extractAfter(func2str(g),')'),a,err(i),37);
    
end

if(choice=='d')
    str_fprime=input("Enter f'(x): ",'s');
    str_fprime=append('@(x)',str_fprime);
    fprime=str2func(str_fprime);
    a = input('Enter initial guess: ');
    err=zeros(1,max_it);
    iter=zeros(1,max_it);
    for i=1:max_it
        a_prev=a;
        a=a-(f(a)/fprime(a));
        err(i) = abs(a-a_prev)*100/abs(a);
        iter(i)=i;
        if(err(i)<max_err)
            break;
        end
    end
    line(iter(1:i),err(1:i),'LineWidth',1.75);
    fplot(f);
    fprintf('The root of %s = 0 is %.10f +- %.10f %c\n',extractAfter(func2str(f),')'),a,err(i),37);
    
end

if(choice=='e')
    a=input('Enter smaller value: ');
    b=input("Enter larger value; ");
    err=zeros(1,max_it);
    iter=zeros(1,max_it);
    for i=1:max_it
        c = b-((f(b)*(b-a))/(f(b)-f(a)));
        err(i)=abs((c-b)/c)*100;
        iter(i)=i;
        if(err(i)<max_err)
            break;
        end
        a=b;
        b=c;
    end
    line(iter(1:i),err(1:i),'LineWidth',1.75);
    fplot(f);
    fprintf('The root of %s = 0 is %.10f +- %.10f %c\n',extractAfter(func2str(f),')'),c,err(i),37);
    
end


