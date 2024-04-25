%
% 
%
function ck = getCk(k, q, dq, c, A0)
ck=eye(3,3);
if k==1
    ck(1:2,3)=-[0 1;-1 0]*(0               -A0{k}*c{k,k});
else
    ck(1:2,3)=-[0 1;-1 0]*(A0{k-1}*c{k-1,k}-A0{k}*c{k,k});
end