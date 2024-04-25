%
% 
%
function sk = getSk(k, q, dq, c, A0)    
    sk=[0;0;1];
    sk(1:2)=[0 1;-1 0]*A0{k}*c{k,k};
