% version from Zengyang
function [xhk, pf] = particle_filter_val(sys, yk, pf, resampling_strategy,param, a)

k = pf.k;
if k == 1
    error('error: k must be an integer greater or equal than 2');
end

%% Initialize variables
% Ns = pf.Ns;                              % number of particles
Ns = 1;
nx = size(pf.particles,1);               % number of states

wkm1 = pf.w(:, k-1);                     % weights of last iteration
if k == 2
    for i = 1:Ns                          % simulate initial particles
        pf.particles(:,i,1) = pf.gen_x0(yk,1,param); % at time k=1 (obs,P,param)
    end
    wkm1 = repmat(1/Ns, Ns, 1);           % all particles have the same weight
end


%% Separate memory
xkm1 = pf.particles(:,:,k-1); % extract particles from last iteration;
xk   = zeros(size(xkm1));     % = zeros(nx,Ns);
wk   = zeros(size(wkm1));     % = zeros(Ns,1);

%% Algorithm 3 of Ref [1]
for i = 1:Ns
    xk(:,i) = sys(k, xkm1(:,i), param, pf.gen_sys_noise());
%     xk(:,i) = sys(k, a, param, pf.gen_sys_noise());
    wk(i) = wkm1(i) * pf.p_yk_given_xk(k, yk, xk(:,i),param);
end

%% Normalize weight vector
if any(isnan(wk))
    disp('weighting nan')
end
wk = wk./sum(wk);

%% Calculate effective sample size: eq 48, Ref 1
Neff = 1/sum(wk.^2);

%% Resampling
% remove this condition and sample on each iteration:
% [xk, wk] = resample(xk, wk, resampling_strategy);
%if you want to implement the bootstrap particle filter
resample_percentaje = 0.50;
Nt = resample_percentaje*Ns;
if Neff < Nt
    disp(['Resampling ...',num2str(k)])
    [xk, wk] = resample(xk, wk, resampling_strategy);
    % {xk, wk} is an approximate discrete representation of p(x_k | y_{1:k})
end

%% Compute estimated state
xhk = zeros(nx,1);
for i = 1:Ns
    xhk = xhk + wk(i)*xk(:,i);
end

%% Store new weights and particles
pf.w(:,k) = wk;
pf.particles(:,:,k) = xk;

return; % bye, bye!!!

%% Resampling function
function [xk, wk, idx] = resample(xk, wk, resampling_strategy)

Ns = length(wk);  % Ns = number of particles

switch resampling_strategy
    case 'multinomial_resampling'
        with_replacement = true;
        idx = randsample(1:Ns, Ns, with_replacement, wk);
    case 'systematic_resampling'
        % this is performing latin hypercube sampling on wk
        edges = min([0 cumsum(wk)'],1); % protect against accumulated round-off
        edges(end) = 1;                 % get the upper edge exact
        u1 = rand/Ns;
        % this works like the inverse of the empirical distribution and returns
        % the interval where the sample is to be found
        [~, idx] = histc(u1:1/Ns:1, edges);
    otherwise
        error('Resampling strategy not implemented')
end

xk = xk(:,idx);                    % extract new particles
wk = repmat(1/Ns, 1, Ns);          % now all particles have the same weight

return;  % bye, bye!!!
