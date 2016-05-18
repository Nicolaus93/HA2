%load data
clear all;
atl_data=load('atlantic.mat');
data=atl_data.atlantic;

%% %%Fit Gumbel Distribution to data

[Beta_hat,Mu_hat]=est_gumbel(data);
T=3*14*100;
rtrn_hat=inv_gumbel(1-1/T,Beta_hat,Mu_hat);

%% Bootstrap
n=500;
B=10000;
bootMu=zeros(1,B);
bootBeta=zeros(1,B);
bootReturn=zeros(1,B);
for b=1:B
    
u=rand(1,n);
y_boot=inv_gumbel(u,Beta_hat,Mu_hat);
%estimate parameters given new data
[bootBeta(b),bootMu(b)] =est_gumbel(y_boot);
%estimate of the returin with new estimate parameters

bootRtrn(b)=inv_gumbel((1-1/T),bootBeta(b),bootMu(b));
end

deltaBeta=sort(bootBeta-Beta_hat);
deltaMu=sort(bootMu-Mu_hat);
deltaRtrn=sort(bootRtrn-rtrn_hat);
alpha=0.05;
%confidence intervals
L_Beta = Beta_hat - deltaBeta(ceil((1 - alpha/2)*B));
U_Beta = Beta_hat - deltaBeta(ceil(alpha*B/2));
L_Mu = Mu_hat - deltaMu(ceil((1 - alpha/2)*B));
U_Mu = Mu_hat - deltaMu(ceil(alpha*B/2)); 
L_rtrn = rtrn_hat - deltaRtrn(ceil((1 - alpha)*B));
U_rtrn = rtrn_hat - deltaRtrn(ceil((alpha)*B));


figure
subplot(3,1,1)       
histogram(deltaBeta)
title('Delta Beta')

subplot(3,1,2)       
histogram(deltaMu)
title('Delta Mu')

subplot(3,1,3) 
histogram(deltaRtrn)
title('Delta 100 -Year Return Value')