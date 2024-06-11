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

PT_dBm = 20;
PT = 10^(PT_dBm/10);
alpha_comm_dB = 0*ones(num_users,1); % w.r.t PT
SINR_constraint_dB = 5*ones(num_users,1);
noise_sensing_dB = -10*ones(num_targets,1); % w.r.t PT
noise_comm_dB = -10*ones(num_users,1); % w.r.t PT
alpha_sensing_dB = -10*ones(num_targets,1); % w.r.t PT

std_noise_comm = 10^(noise_comm_dB./20)*sqrt(PT);
std_noise_sensing = 10^(noise_sensing_dB./20)*sqrt(PT);
P_comm = 10^(SINR_constraint_dB/10).*(std_noise_comm.^2);

%% Comm. channel and target response matrix
comm_loc = (rand(num_users,1)-0.5)*2;                                          % communication user location
alpha_H = sqrt(PT*10^(alpha_comm_dB/10))*exp(1i*2*pi*rand());
H = (alpha_H/sqrt(N_TX))*exp(1i*pi*loc_tx'*comm_loc);              % communication channel

radar_loc = (rand(num_targets,1)-0.5)*2;                                          % radar target location
alpha_targets = sqrt(PT*10^(alpha_sensing_dB/10)).*exp(1i*2*pi*rand(num_targets,1));
G = zeros(N_RX,N_TX);
for i=1:num_targets
    G = G + alpha_targets(i)*exp(1i*pi*loc_rx'*radar_loc(i))*exp(-1i*pi*loc_tx*radar_loc(i));
end

%% Plotting
f_fine = linspace(-1,1,1000);
A_fine = exp(1i*pi*loc_rx'*f_fine);

%% Successive elimination
Tmax = 20;
delta = 2/N_TX;
f_start = -1+(delta/2);
f_end = 1-(delta/2);
confidence=0.95;
t=1;
while true
    f_mid = (f_start+f_end)/2;
    txbeam_1 = get_beamformer(N_TX,PT,f_start,f_mid-(delta/2),H,P_comm);
    rxbeam_1 = get_beamformer(N_TX,PT,f_start,f_mid-(delta/2));

    txbeam_2 = get_beamformer(N_TX,PT,f_mid+(delta/2),f_end,H,P_comm);
    rxbeam_2 = get_beamformer(N_TX,PT,f_mid+(delta/2),f_end);
    
    arm1 = BanditArm(txbeam_1,rxbeam_1,confidence,std_noise_sensing); % need to initialise other parameters
    arm2 = BanditArm(txbeam_2,rxbeam_2,confidence,std_noise_sensing); % need to initialise other parameters
    
    plot(f_fine,abs(txbeam_1'*A_fine),'--b','LineWidth',1.5)
    hold on
    plot(f_fine,abs(rxbeam_1'*A_fine),'-b','LineWidth',1.5)
    plot(f_fine,abs(txbeam_2'*A_fine),'--r','LineWidth',1.5)
    plot(f_fine,abs(rxbeam_2'*A_fine),'-r','LineWidth',1.5)
    xline(radar_loc,'LineWidth',1.5);
    xline(comm_loc,'g','LineWidth',1.5);
    hold off

    max_lcb = 0; 
    min_ucb = Inf;
    keep_arm = 0;
    while max_lcb <= min_ucb
        noise = (std_noise_sensing/sqrt(2))*(randn(N_RX,1) + 1i*randn(N_RX,1));
        arm1 = arm1.pull_arm(G,noise,t);
        t = t+1;

        noise = (std_noise_sensing/sqrt(2))*(randn(N_RX,1) + 1i*randn(N_RX,1));
        arm2 = arm2.pull_arm(G, noise, t);
        t = t+1;

        [max_lcb,keep_arm] = max([arm1.LCB(), arm2.LCB()]);
        min_ucb = min([arm1.UCB(), arm2.UCB()]);
    end
   
    if keep_arm == 1
        f_end = f_mid+(delta/2); 
    else
        f_start = f_mid-(delta/2);
    end

    if t>Tmax
       break
    end
    if f_end-f_start<delta
       break
    end
end

