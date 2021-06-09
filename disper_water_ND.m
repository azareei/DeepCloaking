function [k0,kI] =  disper_water_ND(h,n,alpha)

syms x 
k0 = double(abs(vpasolve(alpha./tanh(x*h) == x,x)));

kI = zeros(1,n);

syms X
counter = 1;
for ii = 1:n
     kI(counter) = vpasolve(alpha == -X*tan(X*h),X,1.1*pi/(2*h)+(ii-1)*pi/h);
    
    counter = counter+1;
end

