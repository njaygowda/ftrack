function [aki] = tvlp_l2(x,p,q)

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
%     variable x_l2(p*(q+1));
%     minimize ( norm(A*x_l2-b,2) );
%     subject to
%         norm(x_l2,1) <= 1;
% cvx_end

x_l2 = Ypu\Yn;
%%% l2 norm
aki = reshape(x_l2,q+1,p);

return;