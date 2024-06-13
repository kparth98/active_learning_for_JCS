function w_norm = get_beamformer(N_ant, PT, f_start, f_end, transition_gap, H_comm, P_comm)
    loc_tx = 0:N_ant-1;
    
    pass_gap = f_end-f_start;
    pass_start = -1 + transition_gap;
    f2 = 1;
    % disp([f1,f1+f_gap,f1+f_gap+delta_f,f2])
    w_base = cfirpm(N_ant-1,[pass_start,pass_start+pass_gap,pass_start+pass_gap+transition_gap,1],[1,1,0,0]).';

    w = w_base.*exp(1j*pi*loc_tx'*(f_start-pass_start));
    w_norm = w*(sqrt(PT)/norm(w));

    if nargin>5
        if abs(w_norm'*H_comm)^2 < P_comm
            h_hat = H_comm/norm(H_comm);
            a1 = h_hat'*w_norm;

            w_perp = w_norm - h_hat*a1;
            w_perp = w_perp./norm(w_perp);
            a2 = w_perp'*w_norm;

            b1_mag = sqrt(PT - (P_comm/(norm(H_comm)^2)));
            b1_phase = angle(a2);
            b1 = b1_mag*exp(1i*b1_phase);

            b2_mag = sqrt(P_comm)/norm(H_comm);
            b2_phase = angle(a1);
            b2 = b2_mag*exp(1i*b2_phase);

            % x2_phase = angle(H_comm'*w_perp) - angle(H_comm'*au_hat);
            % x2 = x2_mag*exp(1i*x2_phase);

            w_norm = b1*w_perp + b2*h_hat;
        end
    end

end