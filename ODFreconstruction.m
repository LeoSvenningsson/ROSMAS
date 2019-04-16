clear all

P2=0.46; % Fill in your P_2 value here
%P4=0.28; % Used for the extra functions. Uncomment them if needed

res=401;
thetaV=linspace(0,pi,res);

%%%%%%%%%%%%%%%%%%%
%calculate Wrapped Lorentzian distribution
%%%%%%%%%%%%%%%%%%%
f = @(gamma) legendreWL_P2(gamma,P2);

initWL= 1;

[gamma,fvalWL] = fminsearch(f,initWL);
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
%calculate most probable distribution
%%%%%%%%%%%%%%%%%%%
% f = @(lambda) legendreMP(lambda,P2,P4);
% 
% initMP= [1.1,4];
% 
% [lambda1and2,fvalMP] = fminsearch(f,initMP);
% MP=exp(lambda1and2(1)*(3*cos(thetaV).^2 - 1)/2 + lambda1and2(2)*(35*cos(thetaV).^4 - 30*cos(thetaV).^2 + 3)/8);
% MP=MP/(integral(@(theta) exp(lambda1and2(1)*(3*cos(theta).^2 - 1)/2 + lambda1and2(2)*(35*cos(theta).^4 - 30*cos(theta).^2 + 3)/8),0,pi));
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
%calculate gauss distribution
%%%%%%%%%%%%%%%%%%%
% f = @(mphi) legendreGauss(mphi,P2,P4);
% 
% initGauss= [0.5,0];
% con1a=[];
% con1b=[];
% con2a=[];
% con2b=[];
% lb=[0,0];
% ub=[10^3,pi/2];
% nonlcon=[];
% options = optimoptions('fmincon','display','none');
% [mandphi,fvalGauss] = fmincon(f,initGauss,con1a,con1b,con2a,con2b,lb,ub,nonlcon,options);
% G=sqrt(mandphi(1)/pi)*exp(-mandphi(1)*(thetaV-mandphi(2)).^2);
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
%plot wrapped lorentzian distribution
%%%%%%%%%%%%%%%%%%%
WL=1/(pi)*sinh(gamma)./(cosh(gamma)-cos(2*thetaV))/trapz(thetaV,1/(pi)*sinh(gamma)./(cosh(gamma)-cos(2*thetaV)));
IntegralWL=trapz(thetaV,WL);
figure(1)
plot(thetaV(1:(res-1)/2 +1)*180/pi,WL(1:(res-1)/2 +1),'LineWidth',4);
xlabel('Angle \theta')
ylabel('ODF f(\theta)')
title('Wrapped Lorentzian ODF')
%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%
%plot all
%%%%%%%%%%%%%%%%%%%
% figure(4)
% plot(thetaV(1:(res-1)/2 +1)*180/pi,WL(1:(res-1)/2 +1),thetaV(1:(res-1)/2 +1)*180/pi,MP(1:(res-1)/2 +1),thetaV(1:(res-1)/2 +1)*180/pi,G(1:(res-1)/2 +1),'LineWidth',4)
% xlabel('Angle \theta')
% ylabel('ODF f(\theta)')
% legend('Wrapped Lorentzian','Most Probable','Gauss')
% title('All ODFs')
%%%%%%%%%%%%%%%%%%%
