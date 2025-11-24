% Material balance based on best model candidate from incremental
% identification

function dc = rate_model_sim(t,c,N,par)

[m,n] = size(N);
dc = zeros(n,1);

r = [par(1)*c(1)*c(2);par(2)*c(2)*c(3);par(3)*c(2)*c(5);]; % reaction rate vector

dc = N'*r; % Material balance equation