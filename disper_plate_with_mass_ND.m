function [mu0, muI, mu1,mu2,mu3,mu4]=disper_plate_with_mass_ND(h,n, alpha, beta, gamma);

syms x 
mu0 = double(abs(vpasolve(alpha./tanh(x*h) == (1-alpha*gamma+beta*x^4)*x,x)));


muI = zeros(1,n);

syms X
counter = 1;

%keyboard
    
for ii = 1:1:n    
    muI(counter) = vpasolve(-alpha./((1-alpha*gamma+beta*X^4))==X*tan(X*h),X,1.1*pi/(2*h)+(ii-1)*pi/h);

    counter = counter+1;
end

tt=(1-alpha*gamma)/beta;
tt2=(sqrt(2*sqrt(tt))+1i*sqrt(2*sqrt(tt)))/2;
tt3=(-sqrt(2*sqrt(tt))+1i*sqrt(2*sqrt(tt)))/2;

mu1 = vpasolve(-alpha./((1-alpha*gamma+beta*X^4))==X*tan(X*h),X,tt3+pi/h*(1-1i)/100);
mu2 = conj(mu1);

mu3 = vpasolve(-alpha./((1-alpha*gamma+beta*X^4))==X*tan(X*h),X,tt2+pi/h*(1-1i)/100);
mu4 = conj(mu3);


