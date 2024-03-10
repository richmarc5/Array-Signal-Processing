
clear all;
close all;
h = [0.5; 1; -0.6]; % Channel to estimate
 mu = 0.01; % Stepsize
 trainingSamples = 1000;
 x = sign(randn(trainingSamples,1)); % Generate BPSK data
 r = filter(h,1,x); % Apply channel
 L = length(h); h_hat = zeros(L,1);
 err=[];
 %% Estimate channel
 for n = L:trainingSamples
 % Select part of training input
 in = x(n:-1:n-L+1);
 % Apply channel estimate to training data
 y = h_hat'*in;
 % Compute error
 e = r(n)-y;
 % Update taps
 h_hat = h_hat + mu*conj(e)*in
%  fprintf('h_hat= %f\n',h_hat)
 fprintf('err= %f\n',e)
 err=[err e];
 end
h
figure(1)
plot(err,'-r')