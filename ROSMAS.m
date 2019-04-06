%Version 1.0
%Licenced by GPLv3
%Free to use share and adapt
%Appropriate credits given to Leo Svenningsson and relevant article(s) with citation

clear all
tic
%%%%%
randN = 1000000;                  % number of Monte Carlo samples for a 5 dimensional integral
err_est = 1;                      % error estimation (0 disabled) (1 enabled) may increase computaion time/memory
%%%%%
L=2;                            % Legendre number
M=1;                            % Rotation slice in the M row
N=1;                            % Number of N peaks calculatd from -N to N
%%%%%
omegaL = 125.8;                 % Larmor frequency MHz
omegaR = 1500;                   % Spinning frequency Hz
%%%%%
beta1 = 60*pi/180;               % Director frame to rotor frame (radians)
%%%%%
sigmaPAS = zeros(3,3);            % Empty CSA tensor matrix
sigmaPAS(1,1) = -122.21;          % Chemical shielding/shift anisotropy sigma_11 (ppm)
sigmaPAS(2,2) = -108.42;          % Chemical shielding/shift anisotropy sigma_22 (ppm)
sigmaPAS(3,3) = -90.96;           % Chemical shielding/shift anisotropy sigma_33 (ppm)
%%%%%
alpha0 = 0*pi/180;                % Euler rotation angle alpha0 (radians) for Principal Axis System to Molecular Frame
beta0 = 0*pi/180;                 % Euler rotation angle beta0 (radians) for Principal Axis System to Molecular Frame
gamma0 = 0*pi/180;                % Euler rotation angle gamma0 (radians) for Principal Axis System to Molecular Frame
%%%%%

%Input ends here

%%%%%
PAStoMF=[alpha0 beta0 gamma0]; % Euler rotation angles (radians) for Principal Axis System to Molecular Frame [alpha_0 beta_0 gamma_0]
%%%%%

rotz1=[cos(PAStoMF(1)) -sin(PAStoMF(1)) 0; sin(PAStoMF(1)) cos(PAStoMF(1)) 0; 0 0 1]; % Z rotation
rotx2=[1 0 0; 0 cos(PAStoMF(2)) -sin(PAStoMF(2)); 0 sin(PAStoMF(2)) cos(PAStoMF(2))]; % X rotation
rotz3=[cos(PAStoMF(3)) -sin(PAStoMF(3)) 0; sin(PAStoMF(3)) cos(PAStoMF(3)) 0; 0 0 1]; % Z rotation

R=rotz1*rotx2*rotz3;            % pre-multiplied rotation matrix

sigmaMF=R*sigmaPAS*transpose(R); % Euler rotation for a second rank tensor

%%%%%
%mu=omegaL*(sigma33-sigma11)/omegaR; % variable, only investigated in
%Harbison et al. 1986
%%%%%

Nlimit=N;
lengthN=length(-Nlimit:Nlimit);
I_LMN=zeros(lengthN,1);
err=zeros(lengthN,1);

for i = 1:2*Nlimit+1
    N=i-Nlimit-1
    [I_LMN(i),err(i)] = calcI(L,M,N,omegaL,omegaR,beta1,sigmaMF,randN,err_est);
end


sumI_MN = sum(I_LMN);

time=toc;
disp(['your calculation took ',num2str(round(time,4,'significant')),' seconds'])

if err_est==0
    disp(['real I_LMN is calculated to:'])
    I_LMN
else
    disp(['I_LMN is calculated to:'])
    I_LMN
    disp(['with a real standard deviation of:'])
    err
end
