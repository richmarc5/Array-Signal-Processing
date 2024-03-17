clear all
close all

clear all
close all
%pkg load signal
% make a normal signal
f1=8e6;
f2=7.5e6;
fs=f1*2.5;
dt=1/fs;
Pw=15e-6;
t=[0:dt:Pw];
y=sin(2*pi*f1*t)+sin(2*pi*f1*t);

% make a nonuniform t domain an make signal from that
trand=0;
for ii=1:length(t)
trand=trand+0.1*dt*rand;% nonuniform t
trandsave(ii)=trand;%save nonuniform t
yrand(ii)=sin(2*pi*f1*trand)+sin(2*pi*f1*trand);%signal
end

figure(1)
plot(trandsave,y)%plot signal made
fTrand=abs(fft(yrand));%take fft
figure(2)
plot(fTrand)%plot fft vs sample pts


fac=0.1;%oversample by 0.5*dt
 [y,fsNEW] = resampleSINC(dt,fac,trandsave,Pw,yrand);
   
figure(3)
plot(y)

yFT=abs(fft(y));
freq=[0:length(yFT)-1]*fsNEW*1/(length(yFT));
figure(4)
plot(freq,yFT)



