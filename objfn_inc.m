% Objective function for incremental identification

function sse = objfn_inc(t,C,rate_out,i,j,par)

rate_pred = ratemodel_inc(C,i,j,par);

sse = sum((rate_out(:,i) - rate_pred).^2);
end