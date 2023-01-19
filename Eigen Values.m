file=input('Enter name of .txt file: ','s');
fid=fopen(file);
inp=fscanf(fid,'%f');
n=inp(1);
max_it=inp(end-2);
max_err=inp(end-1);
closest=inp(end);
A=reshape(inp(2:end-3),n,n);

prompt = '1.Power method\n2.Inverse power method\n3.Inverse power method with shift\n4.QR method\nChoose number corresponding to your choice: ';
choice=input(prompt);

if(choice==1)
    it=0;
    z=ones(n,1);
    for i=1:max_it
        z_new=(A*z)/max(A*z);
        e_new=max(A*z);
        if(i~=1)
            err=abs(e_new-e)*100/norm(e_new);
            if(err<max_err)
                break;
            end
        end
        e=e_new;
        z=z_new;
        it=it+1;
    end
    z=(A*z)/norm(A*z);
    power_f=fopen('Power Method.txt','wt');
    fprintf(power_f,'Eigenvalue:\n%f\nEigenvector:\n',e);
    for i=1:n
        fprintf(power_f,'%f\n',z(i));
    end
    fprintf(power_f,'Iterations:\n%d',it+1);
    fprintf('\nOutput saved to file Power Method.txt\n');
end


if(choice==2)
    A_in=inv(A);
    it=0;
    z=ones(n,1);
    for i=1:max_it
        z_new=(A\z)/max(A\z);
        e_new=max(A\z);
        if(i~=1)
            err=abs(e_new-e)*100/norm(e_new);
            if(err<max_err)
                break;
            end
        end
        e=e_new;
        z=z_new;
        it=it+1;
    end
    z=(A\z)/norm(A\z);
    power_f=fopen('Inverse Power Method.txt','wt');
    fprintf(power_f,'Eigenvalue:\n%f\nEigenvector:\n',1/e);
    for i=1:n
        fprintf(power_f,'%f\n',z(i));
    end
    fprintf(power_f,'Iterations:\n%d',it+1);
    fprintf('\nOutput saved to file Inverse Power Method.txt\n');
end



if(choice==3)
    A_shift=A-closest*eye(n,n);
    it=0;
    z=ones(n,1);
    for i=1:max_it
        z_new=(A_shift\z)/max(A_shift\z);
        e_new=max(A_shift\z);
        if(i~=1)
            err=abs((e_new-e)*100/e_new);
            if(err<max_err)
                break;
            end
        end
        e=e_new;
        z=z_new;
        it=it+1;
    end
    z=(A_shift\z)/norm(A_shift\z);
    power_f=fopen('Inverse Power_Shift.txt','wt');
    fprintf(power_f,'Eigenvalue:\n%f\nEigenvector:\n',closest+(1/e));
    for i=1:n
        fprintf(power_f,'%f\n',z(i));
    end
    fprintf(power_f,'Iterations:\n%d',it+1);
    fprintf('\nOutput saved to file Inverse Power_Shift.txt\n');
end



if(choice==4)
    R=qr(A);
    Q=(A/R);
    err_mat=zeros(n,1);
    it=0;
    for i=1:max_it
        A_new=R*Q;
        R=qr(A_new);
        Q=(A_new/R);
        for j=1:n
                err_mat(j)=abs((A_new(j,j)-A(j,j))/A_new(j,j))*100;
        end

        m = max(err_mat);
        it=it+1;
        if(m<max_err)
            break;
        end
        A=A_new;

    end

    qr_f=fopen('QR Method.txt','wt');
    fprintf(qr_f,'Eigenvalues:\n');
    for i=1:n
        fprintf(qr_f,'%f\n',A(i,i));
    end
    fprintf(qr_f,'Iterations: %d\n',it);
    fprintf('\nOutput saved to file QR Method.txt\n');

end