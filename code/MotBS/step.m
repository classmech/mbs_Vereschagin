function y=step(x,x0,y0,x1,y1)

a1=y1-y0;
delta=(x-x0)./(x1-x0);

y=(x<x0).*y0+((x>=x0).*(x<=x1)).*(y0+a1.*delta.*delta.*(3-2*delta))+(x>x1).*y1;    
