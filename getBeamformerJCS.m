function w=getBeamformerJCS(BS_param,Comm_param,f_targets)
% Based on Yonina Eldar's paper on CRB optimisation

    n_targets = length(f_targets);
    if Comm_param.num_users == 1 && n_targets == 1
        a = exp(1i*pi*BS_param.loc_tx'*f_targets);
        P_T = 10^(Comm_param.PTot_dBm/10);
        N_TX = size(BS_param.A_TX,1);
        Gamma = 10^(Comm_param.SINR_constraint_dB/10);
        noise_comm_std = 10^(Comm_param.noise_comm_dBm/10);
        h = Comm_param.H';

        if P_T*abs(h'*a)^2 > N_TX*Gamma*noise_comm_std^2
            w = sqrt(P_T/N_TX)*a;
        else
            u1 = h./norm(h);
            au = (eye(N_TX) - u1*u1')*a;
            au = au./norm(au);

            temp = Gamma*noise_comm_std^2/norm(h)^2;

            x1 = sqrt(temp)*(u1'*a)/abs(u1'*a);
            x2 = sqrt(P_T - temp)*(au'*a)/abs(au'*a);

            w = x1*u1 + x2*au;
        end
    else
        error('not implemented')
    end
end