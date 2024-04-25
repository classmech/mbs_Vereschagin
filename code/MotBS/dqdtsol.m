function dx = dqdtsol(t,x)
%DQDT Summary of this function goes here
%   Detailed explanation goes here

global Model;

n=(size(x,1)-3)/2;

q =x(1:n);
dq=x(n+1:2*n);

dx=[dq;dq;0;0;0];

dx(7)=step(q(1),      0, 1, +0.001, 0);
dx(8)=step(q(2), -0.001, 0,      0, 1);
dx(9)=step(q(3),      0, 1, +0.001, 0);

A0=cell(n,1); 
A0{1}=getA(q(1));
for i=2:n
	A0{i}=A0{i-1}*getA(q(i));    
end

Q=cell(n,1);

M10=-2.0;
M1K=-1.4;
M20=+1.5;
M2K=+0.5;
M30=-0.75;
M3K=-0.1;

M1=M10-(M10-M1K)/(0.5*pi)*(pi*0.5-q(1))-step(x(7),0,0,0.0001,1)*(q(1)*1e4+dq(1)*100);
M2=M20-(M20-M2K)/(pi)*(pi    +q(2))-step(x(8),0,0,0.0001,1)*(q(2)*1e4+dq(2)*100);
M3=M30-(M30-M3K)/(pi)*(pi    -q(3))-step(x(9),0,0,0.0001,1)*(q(3)*1e4+dq(3)*100);

Q{1}=[0;0;M1-M2];
Q{2}=[0;0;M2-M3];
Q{3}=[0;0;   M3];

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

