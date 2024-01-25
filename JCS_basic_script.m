%% Basic Script for Joint Comm. and Sensing
clear all

%% Define Parameters (following Yonina's paper on CRB minimisation)
N_TX = 16;                      % No. of transmit(TX) antennas at base station
N_RX = 16;                      % No. of receive(RX) antennas at base station
BS_param.noise_sensing_dBm = 0;%in dBm       % noise power (in dBm) at base station RX
Comm_param.noise_comm_dBm = 0;%in dBm        % noise power (in dBm) at comm. user RX

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
Comm_param.num_users = 1;                                   % no. of users (assumed to be known)
Comm_param.loc_idx = randi(N_grid,Comm_param.num_users);          % index of gridpoint where users are located
Comm_param.H = BS_param.A_TX(:,Comm_param.loc_idx)';              % Comm. channel matrix (can be replaced by a random matrix)

Target_param.num_targets = 1;        % no. of targets (can be an upper bound)
Target_param.loc_idx = randi(N_grid,Target_param.num_targets);     % index of gridpoint where targets are located 
Target_param.alpha_stddev = 10^(-10/10)*ones(Target_param.num_targets); % SNR of target is -10dB here
Target_param.alpha = Target_param.alpha_stddev.*(randn(Target_param.num_targets) + ...
                                        1i*randn(Target_param.num_targets))/sqrt(2); % the complex reflection coefficients of each target

for i=1:Target_param.num_targets
    G = Target_param.alpha(i)*BS_param.A_RX(:,Target_param.loc_idx(i))*BS_param.A_TX(:,Target_param.loc_idx(i))';
end


Comm_param.PTot_dBm = 30;% in dBm                            % base station TX power constraint (in dBm)
Comm_param.SINR_constraint_dB = 30*ones(Comm_param.num_users,1);        % minimum SINR for each user for communication

Tmax = N_grid;
%% 
f_fine = linspace(-1,1,10000);
A_TX_fine = exp(1i*pi*BS_param.loc_tx'*f_fine);
noise_std = 10^(BS_param.noise_sensing_dBm/10);

for t=1:Tmax
    f_target = f_grid(t); % sampling each grid point 
    w_t = getBeamformerJCS(BS_param, Comm_param, f_target);

    % beampattern
    % plot(f_fine, 10*log10(abs(w_t'*A_TX_fine)),'LineWidth',1.5);
    plot(f_grid, abs(w_t'*BS_param.A_TX),'LineWidth',1.5);
    hold on
    xline(f_grid);
    xline(f_grid(Target_param.loc_idx),'r','LineWidth',1.5)
    xline(f_grid(Comm_param.loc_idx),'b','LineWidth',1.5)
    ylabel('Beam gain (in dB)')
    xlabel('sin(\theta) grid')
    hold off

    ct = 1; % keeping the comm. symbol constant for now
    x_t = w_t*ct;
    
    noise = noise_std*(randn(N_RX,1) + 1i*randn(N_RX,1))/sqrt(2);
    y_t = G*x_t + noise;

    surrogate_fn = abs(x_t'*y_t);
    disp(surrogate_fn);
end


