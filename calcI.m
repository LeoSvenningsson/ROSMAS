function [I,err] = calcI(L,M,N,omegaL,omegaR,beta1,sigma,randN,err_est)

A1 = @(alpha2,beta2) -sqrt(2)/3*omegaL/omegaR*(sin(2*beta2)*((sigma(2,2) + sigma(1,1))/2 - (sigma(2,2)-sigma(1,1))*cos(2*alpha2)/2 + sigma(1,2)*sin(2*alpha2) - sigma(3,3)) + 2*cos(2*beta2)*(sigma(1,3)*cos(alpha2) + sigma(2,3)*sin(alpha2)));

B1 = @(alpha2,beta2) -2*sqrt(2)/3*omegaL/omegaR*(sin(beta2)*(((sigma(2,2)-sigma(1,1))*sin(2*alpha2)/2) + sigma(1,2)*cos(2*alpha2)) + cos(beta2)*(sigma(2,3)*cos(alpha2)-sigma(1,3)*sin(alpha2)));

A2 = @(alpha2,beta2) -omegaL/(6*omegaR)*((sigma(2,2) + sigma(1,1))*(cos(2*beta2) - 1)/4 - (sigma(2,2) - sigma(1,1))*cos(2*alpha2)*(3 + cos(2*beta2))/4 + sigma(1,2)*sin(2*alpha2)*(3 + cos(2*beta2))/2 - sigma(1,3)*sin(2*beta2)*cos(alpha2) - sigma(2,3)*sin(alpha2)*sin(2*beta2) + sigma(3,3)*(1 - cos(2*beta2))/2);

B2 = @(alpha2,beta2) -omegaL/(3*omegaR)*(cos(beta2)*((sigma(2,2)-sigma(1,1))*sin(2*alpha2)/2 + sigma(1,2)*cos(2*alpha2)) + sin(beta2)*(sigma(1,3)*sin(alpha2) - sigma(2,3)*cos(alpha2)));


F = @(alpha2,nu,beta2) exp(1i*(A2(alpha2,beta2).*sin(2*nu) - B2(alpha2,beta2).*cos(2*nu) + A1(alpha2,beta2).*sin(nu) - B1(alpha2,beta2).*cos(nu)));
FN = @(alpha2,nu,beta2) F(alpha2,nu,beta2).*exp(-1i*N*nu);
conjF = @(alpha2,gamma,beta2) exp(-1i*(A2(alpha2,beta2).*sin(2*gamma) - B2(alpha2,beta2).*cos(2*gamma) + A1(alpha2,beta2).*sin(gamma) - B1(alpha2,beta2).*cos(gamma)));
conjFNM = @(alpha2,gamma,beta2) conjF(alpha2,gamma,beta2).*exp(1i*(N-M)*gamma);

KNM = @(alpha2,nu,gamma,beta2,gamma2) FN(alpha2,nu,beta2).*conjFNM(alpha2,gamma,beta2).*exp(1i*M*gamma2)*1/(4*pi^2);

legPoly = LegendreFunc(beta1,L);

INML = @(alpha2,nu,gamma,beta2,gamma2) KNM(alpha2,nu,gamma,beta2,gamma2).*legPoly(beta2,gamma2).*sin(beta2)*(2*L+1)/(8*pi^2);

IMC = 0;
if err_est==0
parfor i=1:randN  
  randInt = rand(1,5).*[2*pi,2*pi,2*pi,pi,2*pi];
  IMC = feval(INML,randInt(1),randInt(2),randInt(3),randInt(4),randInt(5))*(pi*(2*pi)^4) + IMC;    
end
I = IMC/randN;
err=NaN;
else
IMCV=(zeros(1,randN));
parfor i=1:randN  
  randInt = rand(1,5).*[2*pi,2*pi,2*pi,pi,2*pi];
  IMCV(i)=feval(INML,randInt(1),randInt(2),randInt(3),randInt(4),randInt(5))*(pi*(2*pi)^4)
end
I=mean(IMCV);
err=std(real(IMCV))/sqrt(randN);
end

end

