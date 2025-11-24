% Objective function for simultaneous identification

function sse = objfn_sim(par, t, c_0, N, C_meas, idx, k)

    % Solve the ODE
    [t, C_model] = ode15s(@(t, C) rate_model_sim(t, C, N, par), t, c_0);

    % Check if bootstrap parameters are provided
    if nargin > 5 && ~isempty(idx) && ~isempty(k)
        % Bootstrap case
        C_model = C_model(idx{k}, :); % Use specific indices for bootstrap
        diff = C_meas{k} - C_model;
    else
        % Point estimate case
        diff = C_meas - C_model;
    end

    % Compute the sum of squared errors
    s = diff * diff';
    sse = (norm(s, 'fro'))^2;
end

