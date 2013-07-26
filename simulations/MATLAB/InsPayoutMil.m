function [ Mmil ] = InsPayoutMil( T, N, r1, h, F, alpha, rho, beta, P, gamma, g )
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here

rtilde=r1-h;
ftilde=(rtilde/r1)*F;

dt=sqrt(T/N);
dW = dt*randn(1, N+1);

ds=T/N;
s=(0:ds:T);

m0 = (beta*ftilde - P)/(rho-g);
Mmil = zeros(1, N+1);
Mmil(1) = m0;
Mtemp = m0;

[ Xtrue ] = DeerPopMil( T, N, r1, h, F, alpha );

for i=1:N
    Winc = dW(i);
    Mtemp = Mtemp + dt*(rho*Mtemp - beta*Xtrue(i)+P) - Winc*gamma*Xtrue(i) - .5*gamma*alpha*Xtrue(i)*(Winc*Winc-dt);
    Mmil(i+1) = Mtemp;
end

%plot(s, Mmil, 'bx')
%legend('Mmil')