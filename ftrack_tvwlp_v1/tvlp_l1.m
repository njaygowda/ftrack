function [aki] = tvlp_l1(x,p,q)

m=1;
for n=p+1:length(x);
    for i=1:p;
        for j=0:q;
            Yn(m)=x(n);
            Ypu(m,(i-1)*(q+1)+j+1)=(n-p-1)^j*x(n-i);
        end;
    end;
    m=m+1;
end
Yn=Yn(:);

A=Ypu;
b=Yn;

% cvx_begin
%     variable x_l1(p*(q+1));
%     minimize ( norm(A*x_l1-b,1) );
%     subject to
%        norm(x_l1,1) <= 1;
% cvx_end

[m,n]=size(A)
f = [ zeros(n,1); ones(m,1); ones(m,1) ];
Aeq = [ A, -eye(m), +eye(m) ];
lb = [ -Inf(n,1); zeros(m,1); zeros(m,1) ];
opts=optimoptions('linprog','MaxIter',100,'Diagnostics','on','Display','iter');
[xzz,fval,exitflag] = linprog(f,[],[],Aeq,b,lb,[],[],opts);
% x_l1 = xzz(1:n,:)-xzz(n+1:2*n,:);
x_l1 = xzz(1:n,:);

aki = reshape(x_l1,q+1,p);

%%% l2 norm
%a_l2 = reshape(Ypu\Yn,q+1,p);

return;
