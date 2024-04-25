function  wprim = getWprim(k, q, dq, c, A0)

if k==1       
   wprim=[(dq(k)+  0)^2*A0{k}*c{k,k}                         ;0];
else
   wkp=sum(dq(1:k-1)); 
   wprim=[(dq(k)+wkp)^2*A0{k}*c{k,k}-wkp*wkp*A0{k-1}*c{k-1,k};0]; 
end


end

