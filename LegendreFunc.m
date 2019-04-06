function PLconvert = LegendreFunc(beta1,L)

syms x beta2 gamma2
if L == 0
    PLconvert = @(beta2,gamma2) 1;
elseif L == 2 
   PL=1/2*(3*x^2 - 1);
   PLsub = subs(PL,x,sin(beta1)*sin(beta2)*cos(gamma2) + cos(beta1)*cos(beta2));
   PLconvert = matlabFunction(PLsub);
elseif L == 4
    PL=1/8*(35*x^4 - 30*x^2 +3);
    PLsub = subs(PL,x,sin(beta1)*sin(beta2)*cos(gamma2) + cos(beta1)*cos(beta2));
    PLconvert = matlabFunction(PLsub);
else
    f=(x^2-1)^L;
    g=diff(f,x,L);
    PL=g/(2^L*factorial(L));
    PLsub = subs(PL,x,sin(beta1)*sin(beta2)*cos(gamma2) + cos(beta1)*cos(beta2));
    PLconvert = matlabFunction(PLsub);
end

end
