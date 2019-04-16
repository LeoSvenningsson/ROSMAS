function y = legendreWL(gamma,P2)
int1=0;
int2=pi;


fun = @(theta) sinh(gamma)./(cosh(gamma)-cos(2*theta)).*sin(theta);



fun1 = @(theta) fun(theta).*(3*cos(theta).^2 - 1)/2;


y =  (integral(fun1,int1,int2)./integral(fun,int1,int2) - P2).^2;


integral(fun,int1,int2);
end

