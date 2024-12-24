function [A, Z,obj,iter] = ALPC(X, max, k, m, alpha,beta)

num_view = length(X);
for i = 1:num_view
    A{i} = zeros(size(X{i},1), m*k);
end
R = zeros(k, size(X{i},2));
I = eye(k); Y= [];
for ik = 1:k
    Y = [Y repmat(I(:,ik),1,m)];
end
clear I
P = Y;
for iv = 1:num_view
    [Up,~,Sp] = svd(A{iv}*Y', 'econ');
    U{iv} = Up*Sp';
end
clear Up Sp
MAX_ITER = max;
for iter = 1:MAX_ITER
    %--------------------updata Z-------------------%
    tZ1 = 0; tZ2 = 0;
    for iv = 1:num_view
        tZ1 = tZ1 + A{iv}'*A{iv};
        tZ2 = tZ2 + A{iv}'*X{iv};
    end
    Z = (tZ1 + beta*eye(m*k))\(tZ2 + beta*P'*R);
    %--------------------updata A{iv}-------------------%
    for iv = 1:num_view
        A{iv} = (X{iv}*Z' + alpha*U{iv}*P)/(Z*Z' + alpha*eye(m*k));
    end
    %--------------------updata U{iv}-------------------%
    for iv = 1:num_view
        [Up,~,Sp] = svd(A{iv}*P', 'econ');
        U{iv} = Up*Sp';
    end
    clear Up Sp
    %--------------------updata P-------------------%
    UA = 0;
    for iv = 1:num_view
        UA = UA + U{iv}'*A{iv};
    end
    P = (alpha*eye(size(R*R')) + beta*(R*R'))\(UA + beta*R*Z');
    %--------------------updata R-------------------%
    R = P*P'\P*Z;
    %--------------------updata loss----------------%
    N1 = 0;  N2 = 0; normX = 0;
    for iv = 1:num_view
        N1 = norm(X{iv} - A{iv}*Z,'fro').^2 + N1;
        N2 = N2 + alpha*norm(A{iv} - U{iv}*P,'fro').^2;
        normX = normX + norm(X{iv},'fro').^2;
    end
    normX = normX.^2;
    N3 = beta*norm(Z-P'*R,'fro').^2;
    obj(iter) = (N1 + N2 + N3)./normX;
    if (iter>2) && ( abs((obj(iter-1)-obj(iter)))<1e-8)
        iter
        break;
    end
end
