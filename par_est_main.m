%%%% Main file for parameter estimation and model discrimination using
%%%% incremental identification approach
% N_best, X_best are physically feasible stoichiometry and extents  
% C - noisy concentration data, time - time points of measured concentration

clc
clear all
close all

%% loading data (results from piMF)

load("piMF_res.mat","time","C","N_best","X_best")

X_best(1,:) = 0;
c0 = C(1,:); % Initial conc
C_out = X_best*N_best+(ones(length(time),1)*C(1,:)); % concentration from computed extents
L = length(C_out);
[nR,nS] = size(N_best);

%% computing the rate of reaction from polynomial fitted extents (for rate-based incremental identification)

% fitting extents to polynomial of order from 2 to 5 and choose best fit
% based on AIC value

for j = 2:5
for i = 1:nR
    p{j-1}(i,:) = polyfit(time,X_best(:,i),j);
    X_fit{j-1}(:,i) = polyval(p{j-1}(i,:),time);
    err_fit(j-1,i) = sum((X_best(:,i)-X_fit{j-1}(:,i)).^2);
    AIC_fit(j-1,i) = (2*j - 2*log(err_fit(j-1,i)));
end
end

[~,order_best_ply] = min(AIC_fit,[],1);
order_best_ply = unique(order_best_ply); 

%%% rate of reaction computed from best polynomial fit of extents
for i = 1:nR
dp(i,:) = polyder(p{order_best_ply}(i,:)); 
rate_est_ex(:,i) = polyval(dp(i,:),time); % computed rate from extents 
end

%% Model discrimination using rate-based incremental identification and simultaneous identification
%%% Defining the number of model candidates and no of parameters in each
%%% kinetic model candidates
% From N_best we see that all reactions are bimolecular. In this work we have used zero order,
% first order (in both species), second order (in both species and overall second order) and reversible second
% order reaction as model candidates. 

mod_cand = [7 7 7]; % no of model candidates for each reaction
np = [ones(1,6) 2; ones(1,6) 2; ones(1,6) 2]'; % no of parameters in each model 
options=optimset('MaxIter',1000,'TolFun',1e-10,'TolX',1e-10,'MaxFunEvals',3000,'Algorithm','interior-point');

x0 = 0.1;

for i = 1:nR
    for j = 1:mod_cand(i)
        
        if(np(j,i)==1)
        [y{i}(j,:),fval{i}(j,:)] = fmincon(@(par)objfn_inc(time,C_out,rate_est_ex,i,j,par),ones(1,np(j,i))*x0,[],[],[],[],ones(1,np(j,i))*1e-2,ones(1,np(j,i))*1e2,[],options);
        else
        [yeq{i},fval{i}(j,:)] = fmincon(@(par)objfn_inc(time,C_out,rate_est_ex,i,j,par),ones(1,np(j,i))*x0,[],[],[],[],ones(1,np(j,i))*1e-2,ones(1,np(j,i))*1e2,[],options);    
        end
        AIC{i}(j,:) = 2*np(j,i)+(L*log(fval{i}(j,:)))+((2*np(j,i)*(np(j,i)+1))/(L-np(j,i)-1)); % AICC criteria (for small data sets)
    end
    [m(i),I(i)] = min(AIC{i}); % Compute candidate with min AIC value
    rate_pred(:,i)= ratemodel_inc(C_out,i,I(i),y{i}(I(i))); % rate predictions for individual reactions based on best model candidate
end

%%% simultaneous identification 
% Initial guess for simultaneous identification: parameter estimate from incremental identification for best candidate
x0_sim = [y{1}(I(1)) y{2}(I(2)) y{3}(I(3))]; 

[y_sim,fval] = fmincon(@(par)objfn_sim(par,time,c0,N_best,C),x0_sim,[],[],[],[],1e-2*ones(1,nR),1e2*ones(1,nR),[],options); % point estimate


%%% Bootstrap the concentration data to find 95% confidence interval

for i = 1:1000
[C_b{i},idx{i}] = datasample(C,size(C,1)); % Generate bootstrap matrices

% estimate for each bootstrap
[yb_sim(i,:),fb_val(i,:)] = fmincon(@(par)objfn_sim(par,time,c0,N_best,C_b,idx,i),y_sim,[],[],[],[],1e-2*ones(1,nR),ones(1,nR),[],options);
end

[yhat,I] = sort(yb_sim);

% removing tails from the bootstrap estimates
% first 25 and last 25 elements are removed
ind = setdiff(1:1000,[0.025*1000:0.975*1000]); 
yhat(ind,:) = [];

ci = [mean(yhat)-2.*std(yhat) ;mean(yhat)+2.*std(yhat)]';
yhat_m = mean(yhat);
[time,cb_sim] = ode15s(@(t,c)rate_model_sim(time,c,N_best,yhat_m),time,c0);
