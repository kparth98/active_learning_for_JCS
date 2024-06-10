N_TX = 16;                      % No. of transmit(TX) antennas at base station
N_RX = 16;                      % No. of receive(RX) antennas at base station (Yonina uses N_RX=20)

num_targets = 1;
num_users = 1;

PT_dBm = 10;
SINR_constraint_dB = 5*ones(num_users,1);
noise_sensing_dB = -10; % w.r.t PT
noise_comm_dB = 0; % w.r.t PT
alpha_comm_dB = 0; % w.r.t PT
alpha_sensing_dB = 0*ones(num_targets,1); % w.r.t PT

% Defining grid in sin(theta) domain
N_grid = N_TX;
delta_f = 2/N_grid;
f_grid = -1 + delta_f/2 + (0:N_grid-1).*delta_f; % uniform grid in interval sin(theta) domain [-1,1)
theta_grid = asin(f_grid);                       % theta_i = sin^(-1)(f_i), i=1,...,N_grid

% Defining array manifolds
loc_tx = 0:N_TX-1;                              % base station TX antenna locations
loc_rx = 0:N_RX-1;                              % base station RX antenna locations
A_TX = exp(1i*pi*loc_tx'*f_grid);                % TX array manifold, A_TX = [a_TX(theta_1) ... a_TX(theta_{N_grid})]
A_RX = exp(1i*pi*loc_rx'*f_grid);                % RX array manifold, A_RX = [a_RX(theta_1) ... a_RX(theta_{N_grid})]

PT = 10^(PT_dBm/10);

% Comm channel
comm_loc = (rand()-0.5)*2;
alpha_H = sqrt(PT*10^(alpha_comm_dB/10))*exp(1i*2*pi*rand());
H = (alpha_H/sqrt(N_TX))*exp(1i*pi*loc_tx'*comm_loc)';

% For plotting
f_fine = linspace(-1,1,10000);
A_TX_fine = exp(1i*pi*loc_tx'*f_fine);

%% Generate codebook
W_codebook = getBeamformerCodebook(N_TX, N_grid, PT);

%% Visualise first beam of each layer
for i=1:4
    plot(f_fine, abs(A_TX_fine'*W_codebook{i}(:,1)))
    hold on
%     plot(f_fine, abs(A_TX_fine'*W_codebook{i}(:,2)))
    xlabel('target location ($\sin(\theta)$)','interpreter','latex')
    ylabel('response, $|\mathbf{a}_{\rm Tx}(\theta)^H\mathbf{w}|$','interpreter','latex')
%     fontsize(16,"points")
end

