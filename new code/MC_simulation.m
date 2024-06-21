%% Script for Joint Comm. and Sensing with successive elimination
clear all
% rng(0)
%% Define Parameters
N_TX = 16;                      % No. of transmit(TX) antennas at base station
N_RX = 16;                      % No. of receive(RX) antennas at base station (Yonina uses N_RX=20)

loc_tx = 0:N_TX-1;              % base station TX antenna locations
loc_rx = 0:N_RX-1;              % base station RX antenna locations

num_targets = 1;
num_users = 1;

PT_dBm = 15;
PT = 10^(PT_dBm/10);
alpha_comm_dB = 0*ones(num_users,1); % w.r.t PT
SINR_constraint_dB = 10*ones(num_users,1);
noise_sensing_dB = 0*ones(num_targets,1); % w.r.t PT
noise_comm_dB = -10*ones(num_users,1); % w.r.t PT
alpha_sensing_dB = -10*ones(num_targets,1); % w.r.t PT

std_noise_comm = 10^(noise_comm_dB./20)*sqrt(PT);
std_noise_sensing = 10^(noise_sensing_dB./20)*sqrt(PT);
P_comm = 10^(SINR_constraint_dB/10).*(std_noise_comm.^2);

%% Comm. channel and target response matrix
comm_loc = (rand(num_users,1)-0.5)*2;                                          % communication user location
alpha_H = sqrt(PT*10^(alpha_comm_dB/10))*exp(1i*2*pi*rand());
H = (alpha_H/sqrt(N_TX))*exp(1i*pi*loc_tx'*comm_loc);              % communication channel

% radar_loc = (rand(num_targets,1)-0.5)*2;                                          % radar target location
% alpha_targets = sqrt(PT*10^(alpha_sensing_dB/10)).*exp(1i*2*pi*rand(num_targets,1)); % exact, not random
% G = zeros(N_RX,N_TX);
% for i=1:num_targets
%     G = G + alpha_targets(i)*exp(1i*pi*loc_rx'*radar_loc(i))*exp(-1i*pi*loc_tx*radar_loc(i));
% end

f_fine = linspace(-1,1,1000);
%%
if false
n_order=16;
transition_gap = 3/n_order;
pass_gap = 0.5;
w_test = get_beamformer(n_order,PT,-1,pass_gap-1,transition_gap);

A_fine = exp(1i*pi*(0:n_order-1)'*f_fine);
plot(f_fine,abs(w_test'*A_fine),'LineWidth',1.5)
hold on
end
% xline([-1+transition_gap,0+transition_gap,0+2*transition_gap],'LineWidth',1)

%% Successive elimination

num_trials = 100;
num_snr = 5;
snr_vals = linspace(15,0,num_snr);
prob_success = zeros(num_snr,1);
avg_sample_comp = zeros(num_snr,1);
confidence=0.99;
Tmax=50;
plotting=true;
delta = 2/N_TX;

delta_tolerance = delta*1.01;

for i_snr = 1:num_snr
    noise_sensing_dB = snr_vals(i_snr);
    std_noise_sensing = sqrt(10^(noise_sensing_dB./10)*PT);

    for trial = 1:num_trials
        radar_loc = (rand(num_targets,1)-0.5)*2;                                          % radar target location
        alpha_targets = sqrt(PT*10^(alpha_sensing_dB/10)).*exp(1i*2*pi*rand(num_targets,1)); % magnitude fixed, only random phase
        G = zeros(N_RX,N_TX);
        for i=1:num_targets
            G = G + alpha_targets(i)*exp(1i*pi*loc_rx'*radar_loc(i))*exp(-1i*pi*loc_tx*radar_loc(i));
        end

        [f_start,f_end, n_samples] = SuccRegionElim(N_TX, N_RX, G, H, PT, ...
                                        P_comm, std_noise_sensing, confidence, delta_tolerance, ...
                                        Tmax, radar_loc, comm_loc, plotting);
        if f_end-f_start <= delta_tolerance && radar_loc>=f_start && radar_loc<=f_end
            prob_success(i_snr) = prob_success(i_snr) + 1;
            avg_sample_comp(i_snr) = avg_sample_comp(i_snr) + n_samples;
        end
        
    end

end

avg_sample_comp_norm = avg_sample_comp./prob_success;
prob_success_norm = prob_success/num_trials;


%%
figure
plot(snr_vals, prob_success_norm,'-o','Color','b','LineWidth',1.5)
ylim([0,1])
ylabel("Prob. of Success")
xlabel("noise power (in dB w.r.t PT)")
grid on

figure
plot(snr_vals, avg_sample_comp_norm,'-o','Color','b','LineWidth',1.5)
ylabel("Avg sample complexity (success)")
xlabel("noise power (in dB w.r.t PT)")
grid on