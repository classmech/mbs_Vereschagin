function dx = dqdt(t,x)
%DQDT Summary of this function goes here
%   Detailed explanation goes here

global Model;

n=size(x,1)/2;

q =x(1:n);
dq=x(n+1:2*n);

dx=[dq;dq];

A0=cell(n,1); 
A0{1}=getA(q(1));
for i=2:n
	A0{i}=A0{i-1}*getA(q(i));    
end

Q=cell(n,1);
for i=1:n
    Q{i}=[0;-Model.mass{1}(1,1)*9.81;0];
end

[Mk,Qk,iUk]=getMkQkiUk(q,dq,Model.mass,Model.c,Q,A0);

w=[0;0;0];
for i=1:n
    Ck=getCk(i, q, dq, Model.c, A0);
    Sk=getSk(i, q, dq, Model.c, A0);
    Wp=getWprim(i, q, dq, Model.c, A0);
    dx(n+i)=iUk{i}*Sk'*(Qk{i}-Mk{i}*(Ck*w+Wp));
    w=Ck*w+Sk*dx(n+i)+Wp;
end

end

