fstr=input('Enter the ODE to be solved f=dy/dt\n','s');
fstr=append("@(t,y)",fstr);
f=str2func(fstr);
t0=input('Enter initial value t0: ');
y0=input('Enter initial value y0: ');
tf=input('Enter final value tf: ');
h=input('Enter step size h: ');
choice = input('Choose option corresponding to choice:\n1.Forward Euler\n2.2nd order Runge Kutta(Midpoint Method)\n3.4th order Runge Kutta\n');

if(choice==1)
    n=(tf-t0)/h;
    y=zeros(n,1);
    y(1)=y0+h*f(t0,y0);
    for i=2:n
        y(i)=y(i-1)+h*f(t0+(i-1)*h,y(i-1));
    end
    
    eul=fopen('Euler.txt','w');
    fprintf(eul,'t,y\n%f %f\n',t0,y0);
    for i=1:n
        fprintf(eul,'%f %f\n',t0+i*h,y(i));
    end
    fprintf('Output saved to file "Euler.txt"\n');

    hold on
    plot(t0,y0,'.','markersize',12);
    for i=1:n
        plot(t0+i*h,y(i),'.','markersize',12);
    end
    hold off

end

if(choice==2)
    n=(tf-t0)/h;
    y=zeros(n,1);
    t=zeros(n,1);
    k1=f(t0,y0);
    k2=f(t0+(h/2),y0+k1*(h/2));
    t(1)=t0+h;
    y(1)=y0+h*k2;
    for i=2:n
        t(i)=t(i-1)+h;
        k1=f(t(i-1),y(i-1));
        k2=f(t(i-1)+h/2,y(i-1)+k1*h/2);
        y(i)=y(i-1)+k2*h;
    end

    rk2=fopen('RK_2ndOrder.txt','w');
    fprintf(rk2,'t,y\n%f %f\n',t0,y0);
    for i=1:n
        fprintf(rk2,'%f %f\n',t(i),y(i));
    end
    fprintf('Output saved to file "RK_2ndOrder.txt"\n');

    hold on
    plot(t0,y0,'.','markersize',12);
    for i=1:n
        plot(t(i),y(i),'.','markersize',12);
    end
    hold off
            
end



if(choice==3)
    n=(tf-t0)/h;
    y=zeros(n,1);
    t=zeros(n,1);
    k1=f(t0,y0);    
    k2=f(t0+h/2,y0+k1*h/2);
    k3=f(t0+h/2,y0+k2*h/2);
    k4=f(t0+h,y0+k3*h);
    y(1)= y0+(h/6)*(k1+2*k2+2*k3+k4);
    t(1)=t0+h;
    for i=2:n
        k1=f(t(i-1),y(i-1));
        k2=f(t(i-1)+(h/2),y(i-1)+k1*(h/2));
        k3=f(t(i-1)+(h/2),y(i-1)+k2*(h/2));
        k4=f(t(i-1)+h,y(i-1)+k3*h);
        y(i)= y(i-1)+(h/6)*(k1+2*k2+2*k3+k4);
        t(i)=t(i-1)+h;
    end
    rk4=fopen('RK_4thOrder.txt','w');
    fprintf(rk4,'t,y\n%f %f\n',t0,y0);
    for i=1:n
        fprintf(rk4,'%f %f\n',t(i),y(i));
    end
    fprintf('Output saved to file "RK_4thOrder.txt"\n');

    hold on
    plot(t0,y0,'.','markersize',12);
    for i=1:n
        plot(t(i),y(i),'.','markersize',12);
    end
    hold off

end
