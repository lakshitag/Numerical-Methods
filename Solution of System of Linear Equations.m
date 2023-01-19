file=input('Enter name of .txt file: ','s');
fid=fopen(file);
n=fscanf(fid,'%d',[1 1]);
mat=readmatrix(file);

prompt='1.Gauss elimination (GE; without pivoting)\n2.GE (with partial pivoting)\n3.LU decomposition by Doolittle method (without pivoting)\n4.LU decomposition by Crout method (without pivoting)\n5.Cholesky decomposition (for symmetric positive definite matrix)\nChoose the number corresponding to your choice: ';

choice=input(prompt);
A=mat(:,1:n); %coefficient matrix
b=mat(:,n+1); 

if(choice==1)
    for i=1:n
        for j=i+1:n
            mat(j,:) = mat(j,:)-(mat(j,i)/mat(i,i))*mat(i,:);
        end
    end
    x=zeros(n,1);
    for i=n:-1:1
        x(i) = (mat(i,n+1)-mat(i,i+1:n)*x(i+1:n))/mat(i,i);
    end
    ge_fid=fopen('Gauss Elimination.txt','w');
    fprintf(ge_fid,'X\n');
    for i=1:n
        fprintf(ge_fid,'%f\n',x(i));
    end
    fprintf('\nOutput saved to file Gauss Elimination.txt\n');
end

if(choice==2)
    for i=1:n
        [mx,mx_in]=max(abs(A(i:n,i)));
        temp=A(mx_in,:);
        A(mx_in,:)=A(i,:);
        A(i,:)=temp;
        for j=i+1:n
            mat(j,:) = mat(j,:)-(mat(j,i)/mat(i,i))*mat(i,:);
        end
    end
    x=zeros(n,1);
    for i=n:-1:1
        x(i) = (mat(i,n+1)-mat(i,i+1:n)*x(i+1:n))/mat(i,i);
    end
    gep_fid=fopen('Gauss Elimination_Pivoting.txt','w');
    fprintf(gep_fid,'X\n');
    for i=1:n
        fprintf(gep_fid,'%f\n',x(i));
    end
    fprintf('\nOutput saved to file Gauss Elimination_Pivoting.txt\n');
end

if(choice==3)
    l=eye(n,n);
    u=zeros(n,n);
    
    u(1,:)=A(1,:);
    l(:,1)=A(:,1)/u(1,1);
    for i=2:n
        for j=i:n
            u(i,j)=A(i,j)-(l(i,1:i-1)*u(1:i-1,j));
        end
        for k=i+1:n
            l(k,i) = (A(k,i)-(l(k,1:i-1)*u(1:i-1,i)))/u(i,i);
        end
    end
    y=zeros(n,1);
    x=zeros(n,1);
    y(1)=b(1);
    for i=2:n
        y(i)=b(i)-l(i,1:i)*y(1:i);
    end
    x(n)=y(n)/u(n,n);
    for i=n-1:-1:1
        x(i)=(y(i)-(u(i,i+1:n)*x(i+1:n)))/u(i,i);
    end

    doo_fid=fopen('Doolittle.txt','wt');
    fprintf(doo_fid,'L\n');
    for i=1:n
        fprintf(doo_fid,'%f ',l(i,:));
        fprintf(doo_fid,'\n');
    end
    fprintf(doo_fid,'U\n');
    for i=1:n
        fprintf(doo_fid,'%f ',u(i,:));
        fprintf(doo_fid,'\n');
    end
    fprintf(doo_fid,'X\n');
    for i=1:n
        fprintf(doo_fid,'%f\n',x(i));
    end
    fprintf('\nOutput saved to file Doolittle.txt\n');

end

if(choice==4)
    u=eye(n,n);
    l=zeros(n,n);
    l(:,1)=A(:,1);
    u(1,:)=A(1,:)/l(1,1);
    for i=2:n
        for j=i:n
            l(j,i)=A(j,i)-l(j,1:i-1)*u(1:i-1,i);
        end
        for k=i+1:n
            u(i,k) = (A(i,k) - l(i,1:i-1)*u(1:i-1,k))/l(i,i);
        end
    end

    y=zeros(n,1);
    x=zeros(n,1);
    y(1)=b(1)/l(1,1);

    for i=2:n
        y(i)=(b(i)-l(i,1:i)*y(1:i))/l(i,i);
    end
    x(n)=y(n);
    for i=n-1:-1:1
        x(i)=y(i)-(u(i,i+1:n)*x(i+1:n));
    end

    crout=fopen('Crout.txt','wt');
    fprintf(crout,'L\n');
    for i=1:n
        fprintf(crout,'%f ',l(i,:));
        fprintf(crout,'\n');
    end
    fprintf(crout,'U\n');
    for i=1:n
        fprintf(crout,'%f ',u(i,:));
        fprintf(crout,'\n');
    end
    fprintf(crout,'X\n');
    for i=1:n
        fprintf(crout,'%f\n',x(i));
    end
    fprintf('\nOutput saved to file Crout.txt\n');

end

if(choice==5)
    l=chol(A)';
    u=l';
    x=zeros(n,1);
    y=zeros(n,1);
    y(1)=b(1)/l(1,1);
    for i=2:n
        y(i)=(b(i)-(l(i,1:i-1)*y(1:i-1)))/l(i,i);
    end
    x(n)=y(n)/u(n,n);
    for i=n-1:-1:1
        x(i)=(y(i)-(u(i,i+1:n)*x(i+1:n)))/u(i,i);
    end

    ch_fid=fopen('Cholesky.txt','w');
    fprintf(ch_fid,'X:\n');

    for i=1:n
        fprintf(ch_fid,'%f\n',x(i));
    end
    fprintf(ch_fid,'L:\n');
    for i=1:n
        fprintf(ch_fid,'%f ',l(i,:));
        fprintf(ch_fid,'\n');
    end
    fprintf('\nOutput saved to file Cholesky.txt\n');
    
end
