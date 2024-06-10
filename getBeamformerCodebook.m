function W_codebook = getBeamformerCodebook(N_TX, N_grid, PT)
    loc_tx = 0:N_TX-1;
    num_levels = ceil(log2(N_grid));
    W_codebook = cell(num_levels,1);
    delta_f = 2/N_grid;
    f_grid = -1 + delta_f/2 + (0:N_grid-1).*delta_f;

    epsilon=0;
    for l=1:num_levels
        W_mat = zeros(N_TX,2^l);
        if l==num_levels
            epsilon = 1e-2;
        end
        w_filter = cfirpm(N_TX-1,[f_grid(1),f_grid(N_grid/2^l)+epsilon,f_grid(1+(N_grid/2^l)),f_grid(end)],[1,1,0,0]);
        w_filter = w_filter.'*(sqrt(PT)/norm(w_filter));
        W_mat(:,1) = w_filter;

        for k=2:2^l
            W_mat(:,k) = w_filter.*exp(1j*pi*loc_tx'*(f_grid(1 + (k-1)*(N_grid/2^l))-f_grid(1)));
        end
        W_codebook{l} = W_mat;
    end

end