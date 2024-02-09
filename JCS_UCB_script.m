%% Basic Script for Joint Comm. and Sensing
clear all

%% Define Parameters (following Yonina's paper on CRB minimisation)
N_TX = 16;                      % No. of transmit(TX) antennas at base station
N_RX = 16;                      % No. of receive(RX) antennas at base station
BS_param.noise_sensing_dBm = 0;%in dBm       % std. dev of noise at base station RX
Comm_param.noise_comm_dBm = 0;%in dBm          % std. dev of noise at comm. user RX

%% Defining grid in sin(theta) domain
N_grid = N_TX;
delta_f = 2/N_grid;
f_grid = -1 + delta_f/2 + (0:N_grid-1).*delta_f; % uniform grid in interval [-1,1)
theta_grid = asin(f_grid);                       % theta_i = sin^(-1)(f_i), i=1,...,N_grid

%% Defining array manifolds
BS_param.loc_tx = 0:N_TX-1;                              % base station TX antenna locations
BS_param.loc_rx = 0:N_RX-1;                              % base station RX antenna locations
BS_param.A_TX = exp(1i*pi*BS_param.loc_tx'*f_grid);                % TX array manifold, A_TX = [a_TX(theta_1) ... a_TX(theta_{N_grid})]
BS_param.A_RX = exp(1i*pi*BS_param.loc_rx'*f_grid);                % RX array manifold, A_RX = [a_RX(theta_1) ... a_RX(theta_{N_grid})]

%% Setting up scenario and comm. constraints
Comm_param.num_users = 1;                      % no. of users (assumed to be known)
user_loc_idx = randi(N_grid,Comm_param.num_users);         % index of gridpoint where users are located
Comm_param.H = BS_param.A_TX(:,user_loc_idx)';          % Comm. channel matrix (can be replaced by a random matrix)

Target_param.num_targets = 1;        % no. of targets (can be an upper bound)
Target_param.loc_idx = randi(N_grid,Target_param.num_targets);     % index of gridpoint where targets are located 
Target_param.alpha_stddev = ones(Target_param.num_targets);                 % 
Target_param.alpha = Target_param.alpha_stddev.*(randn(Target_param.num_targets) + ...
                                        1i*randn(Target_param.num_targets))/sqrt(2); % the complex reflection coefficients of each target

for i=1:Target_param.num_targets
    G = Target_param.alpha(i)*BS_param.A_RX(:,Target_param.loc_idx(i))*BS_param.A_TX(:,Target_param.loc_idx(i))';
end


Comm_param.PTot_dBm = 30;% in dBm                            % base station TX power constraint (in dBm)
Comm_param.SINR_constraint_dB = 30*ones(Comm_param.num_users,1);        % minimum SINR for each user for communication

Tmax = 50;
%% 
f_fine = linspace(-1,1,10000);
A_TX_fine = exp(1i*pi*BS_param.loc_tx'*f_fine);
noise_std = 10^(BS_param.noise_sensing_dBm/10);

%% Define UCB algorithm params:
beta = 3;
delta = 0.05;
UCB_grid = 1e6 * ones(N_grid, 1);
initial_guess_theta_idx = randi(N_grid);
real_theta_idx = Target_param.loc_idx;
visit_counts = zeros(N_grid, 1);
avg_rewards = zeros(N_grid, 1);

%% Compute optimal surrogate
opt_wt = getBeamformerJCS(BS_param, Comm_param, f_grid(real_theta_idx));
ct = 1; % keeping the comm. symbol constant for now
x_t = opt_wt*ct;
y_t = G*x_t;
opt_surrogate = abs(x_t'*y_t);
disp(opt_surrogate)

%% UCB algorithm
cumu_regret_array = zeros(Tmax, 1);

for t=1:Tmax
    [idx, UCB_grid] = get_ucb(avg_rewards, visit_counts, beta, delta);
    f_target = f_grid(idx); % sampling each grid point 
    w_t = getBeamformerJCS(BS_param, Comm_param, f_target);

    ct = 1; % keeping the comm. symbol constant for now
    x_t = w_t*ct;
    
    noise = noise_std*(randn(N_RX,1) + 1i*randn(N_RX,1))/sqrt(2);
    y_t = G*x_t + noise;

    surrogate_fn = abs(x_t'*y_t);
    if t == 1
        cumu_regret_array(t) = opt_surrogate - surrogate_fn;
    else
        cumu_regret_array(t) = cumu_regret_array(t-1) + opt_surrogate - surrogate_fn;
    end
    
    % update counts
    visit_count_current_idx = visit_counts(idx);
    visit_counts(idx) = visit_count_current_idx + 1;
    avg_rewards(idx) = (avg_rewards(idx) * visit_count_current_idx ...
        + surrogate_fn) / (visit_count_current_idx + 1);
end

plot(cumu_regret_array)


