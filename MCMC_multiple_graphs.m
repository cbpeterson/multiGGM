function [C_save, Sig_save, adj_save, Theta_save, ar_gamma, ar_theta, ...
    nu_save, ar_nu] = MCMC_multiple_graphs(Theta, b_prior, D_prior, n, S, ...
    C, nu, alpha, beta, a, b, my_w, burnin, nmc, disp)
% Modified from code provided by Hao Wang to allow inference on multiple graphs
% Input parameters:
%   Theta: initial value for graph similarity matrix
%   b_prior: vector of degrees of freedom parameters for G-Wishart priors for each group
%   D_prior: p x p x K array of location parameters for G-Wishart prior for each group
%   n: vector of sample sizes for each group
%   S: p x p x K array of sample covariance matrix for each group
%   C: p x p x K array of initial precision matrices for each group
%   nu: p x p matrix of initial values for nu
%   alpha: shape parameter for Theta slab
%   beta: rate parameter for Theta slab
%   a: first parameter for prior on nu
%   b: second parameter for prior on nu
%   my_w: Bernoulli prior parameter for graph similarity indicators
%   burnin: number of burn-in iterations
%   nmc: number of MCMC sample iterations saved
%   disp: T/F for whether to print progress to screen
% Output parameters:
%   C_save: p x p x K x nmc sample of precision matrix
%   Sig_save: p x p x K x nmc sample of covariance matrix
%   adj_save: p x p x K x nmc sample of adjacency matrix
%   Theta_save: K x K x nmc sample of graph similarity matrix
%   ar_gamma: acceptance rate for between-model moves
%   ar_theta: acceptance rate for within-model moves
%   nu_save: p x p x nmc sample of edge-specific parameters
%   ar_nu: acceptance rate for nu

% K is number of sample groups
K = size(Theta, 1);

% p is number of variables
[p] = size(D_prior, 1);

% Posterior parameters of G-Wishart for each group
b_post = b_prior + n;
D_post = D_prior + S;

% Set up matrices for return values
C_save = zeros(p, p, K, nmc);
Sig_save = C_save;
adj_save = C_save;
Theta_save = zeros(K, K, nmc);
nu_save = zeros(p, p, nmc);

% Acceptance rate for gamma based on number of between model moves
ar_gamma = zeros(K, K);

% Acceptance rate for theta based on number of within model moves
n_within_model = zeros(K, K);
ar_theta = zeros(K, K);
ar_nu = zeros(p, p);

% Initial adjacency matrices as defined by intial precision matrices
adj = abs(C) > 1e-5;

% Proposal parameters for MH steps
alpha_prop = 1;
beta_prop = 1;
a_prop = 2;
b_prop = 4;

% Elementwise product zeros out elements of C that are close to 0
C = C .* adj;

% Perform MCMC sampling
for iter = 1:burnin + nmc
    if (disp && mod(iter, 500) == 0)
        fprintf('iter = %d\n', iter);
    end

    % Update graph and precision matrix for each group using code from Hao Wang
    for cur_graph = 1:K

        % Sample off-diagonal elements for each graph
        for i = 1:p-1
            for j = i+1:p
                % Step 1(a) following modification at the bottom of p. 186
                % Bernoulli proposal with odds p(G') * H(e, Sigma) / P(G)
                % Log of odds
                w = log_H(b_prior(:, cur_graph), D_prior(:, :, cur_graph), ...
                    n(cur_graph), S(:, :, cur_graph), C(:, :, cur_graph), i, j) - ...
                    (nu(i, j) + 2 * Theta(cur_graph, :) * squeeze(adj(i, j, :)));
                
                % 1 / (exp(log odds) + 1) = 1 - Bernoulli probability
                w = 1 / (exp(w) + 1);
                
                % current_ij indicates whether edge (i,j) is in G
                current_ij = adj(i, j, cur_graph);
                
                % propose_ij indicates whether edge (i,j) is in G'
                % Proposal will be 1 if rand(1) < w i.e. will be an edge
                % with probability w
                propose_ij = rand(1) < w;
                
                if (propose_ij ~= current_ij)
                    % Step 1(b)
                    [C_prop] = GWishart_NOij_Gibbs(b_prior(:, cur_graph), ...
                        D_prior(:, :, cur_graph), adj(:, :, cur_graph), ...
                        C(:, :, cur_graph), i, j, propose_ij, 0, 1);
                    
                    % Step 2(b) from expression at top of p. 189
                    r2 = log_GWishart_NOij_pdf(b_prior(:, cur_graph), ...
                        D_prior(:, :, cur_graph), C_prop, i, j, current_ij) ...
                        - log_GWishart_NOij_pdf(b_prior(:, cur_graph), ...
                        D_prior(:, :, cur_graph), C_prop, i, j, propose_ij);
                    
                    % Acceptance rate alpha = min(1, e^r2)
                    if (log(rand(1)) < r2)
                        adj(i, j, cur_graph) = propose_ij;
                        adj(j, i, cur_graph) = propose_ij;
                        current_ij = propose_ij;
                    end
                end
                
                % Step 2(c)
                [C(:, :, cur_graph)] = GWishart_NOij_Gibbs(b_post(:, cur_graph), ...
                    D_post(:, :, cur_graph), adj(:, :, cur_graph), ...
                    C(:, :, cur_graph), i, j, current_ij, 0, 0);
            end
        end
        
        % Step 3: Update C and Sigma given graph
        [C(:, :, cur_graph), Sig] = ...
            GWishart_BIPS_maximumClique(b_post(:, cur_graph), ...
            D_post(:, :, cur_graph), adj(:, :, cur_graph), C(:, :, cur_graph), 0, 1);

        if iter > burnin
            Sig_save(:, :, cur_graph, iter-burnin) = Sig;
        end
    end

    % Update the parameters for network relatedness
    for k = 1:K-1
        for m = k+1:K
            % Between model move
            if Theta(k, m) == 0
                theta_prop = gamrnd(alpha_prop, beta_prop);
            else
                theta_prop = 0;
            end
            
            Theta_prop = Theta;
            Theta_prop(k, m) = theta_prop;
            Theta_prop(m, k) = theta_prop;
            
            % Get terms that are a sum over all edges on log scale
            sum_over_edges = 0;
            for i = 1:p-1
                for j = i+1:p
                    sum_over_edges = sum_over_edges + ...
                        log(calc_mrf_C(Theta, nu(i, j))) + 2 * ...
                        (theta_prop - Theta(m, k)) * adj(i, j, k) * adj(i, j, m) - ...
                        log(calc_mrf_C(Theta_prop, nu(i, j)));
                end
            end
            
            % Calculate MH ratio on log scale
            if theta_prop == 0
                log_ar = alpha_prop * log(beta_prop) - log(gamma(alpha_prop)) + ...
                    log(gamma(alpha)) - alpha * log(beta) - ...
                    (alpha - alpha_prop) * log(Theta(m, k)) + ...
                    (beta - beta_prop) * (Theta(m, k)) + sum_over_edges + ...
                    log(1 - my_w) - log(my_w);
            else
                log_ar = alpha * log(beta) - log(gamma(alpha)) + ...
                    log(gamma(alpha_prop)) - alpha_prop * log(beta_prop) - ...
                    (alpha - alpha_prop) * log(theta_prop) - ...
                    (beta - beta_prop) * theta_prop + sum_over_edges + ...
                    log(my_w) - log(1 - my_w);
            end
            
            % Accept proposal with given probability
            if log_ar > log(unifrnd(0,1))
                Theta(m, k) = theta_prop;
                Theta(k, m) = theta_prop;
                
                % Increment acceptance rate
                ar_gamma(k, m) = ar_gamma(k, m) + 1 / (burnin + nmc);
            end
            
            % Within model model
            if Theta(k, m) ~= 0
                n_within_model(k, m) = n_within_model(k, m) + 1;
                theta_prop = gamrnd(alpha_prop, beta_prop);
                Theta_prop = Theta;
                Theta_prop(k, m) = theta_prop;
                Theta_prop(m, k) = theta_prop;
                
                % Get terms that are a sum over all edges on log scale
                sum_over_edges = 0;
                for i = 1:p-1
                    for j = i+1:p
                        sum_over_edges = sum_over_edges + ...
                            log(calc_mrf_C(Theta, nu(i, j))) + 2 * ...
                            (theta_prop - Theta(m, k)) * adj(i, j, k) * adj(i, j, m) - ...
                            log(calc_mrf_C(Theta_prop, nu(i, j)));
                    end
                end
                
                % Calculate MH ratio on log scale
                log_theta_ar = (alpha - alpha_prop) * (log(theta_prop) - log(Theta(m, k))) + ...
                    (beta - beta_prop) * (Theta(m, k) - theta_prop) + sum_over_edges;
                
                % Accept proposal with given probability
                if log_theta_ar > log(unifrnd(0,1))
                    Theta(m, k) = theta_prop;
                    Theta(k, m) = theta_prop;
                    
                    % Track number of proposals accepted
                    ar_theta(k, m) = ar_theta(k, m) + 1;
                end
            end
        end
    end
    
    % Generate independent proposals for q from beta(a_prop, b_prop) density
    for i = 1:p-1
        for j = i+1:p
            q = betarnd(a_prop, b_prop);
            nu_prop = log(q) - log(1-q);
            
            % Calculate MH ratio on log scale
            % log(p(nu_prop)) - log(p(nu)) + log(q(nu)) - log(q(nu_prop))
            log_nu_ar = (nu_prop - nu(i, j)) * (sum(adj(i, j, :)) + a - a_prop) - ...
                (a + b - a_prop - b_prop) * log(1 + exp(nu_prop)) - ...
                log(calc_mrf_C(Theta, nu_prop)) + ...
                (a + b - a_prop - b_prop) * log(1 + exp(nu(i, j))) + ...
                log(calc_mrf_C(Theta, nu(i, j)));
            
            if log_nu_ar > log(unifrnd(0,1))
                nu(i, j) = nu_prop;
                nu(j, i) = nu_prop;
                ar_nu(i, j) = ar_nu(i, j) + 1 / (burnin + nmc);
            end
        end
    end
    
    % Retain values for posterior sample
    if iter > burnin
        C_save(:, :, :, iter-burnin) = C(:, :, :);
        adj_save(:, :, :, iter-burnin) = adj(:, :, :);
        nu_save(:, :, iter-burnin) = nu;
        Theta_save(:, :, iter-burnin) = Theta;
    end
end

% Compute acceptance rate for theta as number of proposals accepted divided
% by number of within model moves proposed
for k = 1:K-1
    for m = k+1:K
        ar_theta(k, m) = ar_theta(k, m) / n_within_model(k, m);
    end
end