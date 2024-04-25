%
% Структура, содержащая параметры механической системы
%
global Model;
Model=struct;
% Количество тел
n=3;
Model.n=n;
% Шарнирные векторы
% Шарнирный вектор $c_{i,j}$ начинается в центре масс тела и заканчивается
% в шарнире $j$.
Model.c=cell(n,n);
c=cell(n,n);
for i=1:n
    for j=1:n
        c{i,j}=[0;0];        
    end
    c{i,i}=[-0.5-0.005;(-1)^(i+1)*0.02];
    if i<n
        c{i,i+1}=[+0.5+0.005;(-1)^(i)*0.02];
    end
end
Model.c=c;
% Матрицы масс тел системы. Рассматриваемая система плоская, поэтому
% матрица масс каждого тела имеет вид:
% <latex>
% $$ M_i = \begin{pmatrix} m_i & 0 & 0 \\ 0 & m_i & 0 \\ 0 & 0 & J_{iz}
% \end{pmatrix}$$
% </latex> 
% 
Model.mass=cell(n,1);
mass=cell(n,1);
for i=1:n
    mass{i}=[10 0 0; 0 10 0; 0 0 0.5];
end
Model.mass=mass;
% Начальные условия
q0=zeros(n*2+3,1);
for i=1:n
    q0(i)=(-1)^(i+1)*pi;
end
q0(1)=pi/2;
%%
% Интегрирование
[t,q]=ode113(@dqdtsol,[0 8],q0);

%% Запись результатов в файл csv для дальнейшей обработки
dlmwrite('result.csv', [t q(:,1:6)], ',');

%% Построение графиков изменения углов поворота
figure
hold on
plot(t,q(:,1),'k-','LineWidth',2);
plot(t,q(:,2),'k--','LineWidth',2);
plot(t,q(:,3),'k-.','LineWidth',2);
hold off
legend('\phi_1','\phi_2','\phi_3');

%% Построение графиков изменения угловых скоростей
figure
hold on
plot(t,q(:,4),'k-','LineWidth',2);
plot(t,q(:,5),'k--','LineWidth',2);
plot(t,q(:,6),'k-.','LineWidth',2);
hold off
text('interpreter','latex');
legend('\dot \phi_2','\dot \phi_2','\dot \phi_3');
%% Анимация
points=[0.005 1.005 1.005 0.005 0.005; -0.005 -0.005 -0.035 -0.035 -0.005];

M = moviein(size(t,1));

axis([-0.1 3.2 -1.2 2]);

hold on;

p=cell(n,1);

maxY=0;

tmax=max(t);

for i=1:tmax/500:tmax
    cla; 
    a=0;
    dr=[0;0];
    for j=1:3 
        f=interp1(t,q(:,j),i);
        a=a+f;
        
        pp=getA(a)*[points(1,:); points(2,:)*(-1)^(j+1)];
        
        if j>1
            dr=dr+getA(a-f)*[1;-0.04*(-1)^j];
        end
        
        line(pp(1,:) + dr(1), pp(2,:) + dr(2));    
    end
    
    f=getframe;

    % M(:,i) = getframe;    
    % imwrite(M(:,i),'img.png','PNG','Resolution',100);
    % mo= getframe;
    % mov=addframe(mov,mo);               
end

%%
% Energy

qabs=q*0;
v=zeros(size(q,1),2);
h=zeros(size(q,1),2);

T=0;
P=0;

for i=1:n
    v=zeros(size(q,1),2);
    h=zeros(size(q,1),2);
    for j=1:i
       qabs(:,i)=qabs(:,i)+q(:,j); 
       qabs(:,n+i)=qabs(:,n+i)+q(:,n+j);       
    end
    for j=1:i
       if j==i 
           b=-c{j,j};
       else
           b=-c{j,j}+c{j,j+1};
       end
       h=h+reshape(getA(qabs(:,j))*b,size(q,1),2);
       v=v-reshape(getA(qabs(:,j))*[0 1;-1 0]*b,size(q,1),2).*[qabs(:,n+j) qabs(:,n+j)];      
    end    
    T=T+mass{i}(1,1)*(v(:,1).^2+v(:,2).^2)*0.5+mass{i}(3,3)*qabs(:,n+i).^2*0.5;
    P=P+Model.mass{i}(1,1)*9.81*h(:,2);
end

figure; 
E0=T+P;
plot(t,T+P-E0(1));
title('Energy error');




