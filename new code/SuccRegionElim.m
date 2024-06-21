function [f_start,f_end, num_samples] = SuccRegionElim(N_TX, N_RX, G, H, PT, P_comm, std_noise_sensing, confidence, delta_region, Tmax, radar_loc, comm_loc, plotting)
    
    % delta = 2/N_TX;
    % f_start = -1+(1/N_TX);
    % f_end = 1-(delta/2);
    f_start = -1;
    f_end = 1;

    f_fine = linspace(-1,1,1000);
    A_fine = exp(1i*pi*(1:N_TX)'*f_fine);
    
    t=1;
    while true
        f_mid = (f_start+f_end)/2;
        % transition_gap = (f_end-f_start)/N_TX;
        transition_gap = (2)/N_TX;
    
        txbeam_1 = get_beamformer(N_TX,PT,f_start,f_mid,transition_gap,H,P_comm);
        rxbeam_1 = get_beamformer(N_RX,PT,f_start,f_mid,transition_gap);
    
        txbeam_2 = get_beamformer(N_TX,PT,f_mid,f_end,transition_gap,H,P_comm);
        rxbeam_2 = get_beamformer(N_RX,PT,f_mid,f_end,transition_gap);
        
        arm1 = BanditArm(txbeam_1,rxbeam_1,confidence,std_noise_sensing); % need to initialise other parameters
        arm2 = BanditArm(txbeam_2,rxbeam_2,confidence,std_noise_sensing); % need to initialise other parameters
        
        if plotting
        % plot(f_fine,abs(txbeam_1'*A_fine),'--b','LineWidth',1.5)
        % hold on
        plot(f_fine,abs(rxbeam_1'*A_fine),'-','LineWidth',1.5)
        hold on
        ylim([0,20])
        % plot(f_fine,real(rxbeam_1'*A_fine),'--','LineWidth',1.5)
        % hold on
        % plot(f_fine,imag(rxbeam_1'*A_fine),'--','LineWidth',1.5)
        % hold on
        % plot(f_fine,abs(txbeam_2'*A_fine),'--r','LineWidth',1.5)
        plot(f_fine,abs(rxbeam_2'*A_fine),'-','LineWidth',1.5)
        xline(radar_loc,'LineWidth',1.5);
        xline(comm_loc,'g','LineWidth',1.5);
        xregion(f_start,f_end,'FaceColor','r','FaceAlpha',0.15);
        % xline(f_start,'r','LineWidth',1.5);
        hold off
        disp(t)
        end
    
        max_lcb = 0; 
        min_ucb = Inf;
        keep_arm = 0;
        while max_lcb <= min_ucb
            reward_complex_1 = get_measurement(arm1.tx_beam, arm1.rx_beam, G, std_noise_sensing);
            reward_1 = abs(reward_complex_1);
            arm1 = arm1.update_arm(reward_1,t);
            
            if t+1>Tmax
               break
            end
            t = t+1;
    
            reward_complex_2 = get_measurement(arm2.tx_beam, arm2.rx_beam, G, std_noise_sensing);
            reward_2 = abs(reward_complex_2);
            arm2 = arm2.update_arm(reward_2, t);
            
            if t+1>Tmax
               break
            end
            t = t+1;
    
            [max_lcb,keep_arm] = max([arm1.LCB(), arm2.LCB()]);
            min_ucb = min([arm1.UCB(), arm2.UCB()]);
        end
        if t>=Tmax
            break
        elseif f_end-f_start<delta_region
           break
        end
        if keep_arm == 1
            % f_end = f_mid+(transition_gap/2); 
            f_end = f_mid;
        else
            % f_start = f_mid-(transition_gap/2);
            f_start = f_mid;
        end
    end
    num_samples=t;
    if plotting
    disp(num_samples);
    end
end

function y = get_measurement(tx_beam, rx_beam, G, std_noise)
    N_RX = length(rx_beam);
    noise = (std_noise/sqrt(2))*(randn(N_RX,1) + 1i*randn(N_RX,1));
    y = rx_beam'*(G*tx_beam + noise);
end