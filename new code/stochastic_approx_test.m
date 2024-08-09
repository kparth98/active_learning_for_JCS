

N=16;
f_fine = linspace(-4/N,4/N,100);
A_fine = exp(1i*pi*(0:N-1)'*f_fine);

% target_loc = rand()*(2/N)-(1/N);
target_loc = -(1/(N+1));

a_target = exp(1i*pi*(0:N-1)'*target_loc);

shift = 0;
phi1 = 0;
w_test_init = exp(-1i*pi/2)*exp(1i*pi*(0:N-1)'*phi1);
step_size = 1/50;
noise_std = 1;

error_abs = zeros(20,1);

for t=1:20
    w_test = w_test_init.*exp(1i*pi*(0:N-1)'*shift);

    plot(f_fine,real(w_test'*A_fine),'LineWidth',1.5)
    hold on
    grid on
    xline(target_loc,'r','LineWidth',1.5)
    xline(1/N,'LineWidth',1.5)
    xline(-1/N,'LineWidth',1.5)
    hold off

    y = real( w_test'*a_target + noise_std*randn());

    shift_temp = shift - (step_size/t)*y;
    if shift_temp>1/N
        shift = 1/N;
    elseif shift_temp<-1/N
        shift = -1/N;
    else
        shift = shift_temp;
    end

    error_abs(t) = abs(shift - target_loc);
end
% w_test = get_beamformer(N,PT,-1,pass_gap-1,transition_gap);
% A_fine = exp(1i*pi*(0:N-1)'*f_fine);
% phi1 = -2/N;
% w_test1 = exp(1i*pi*(0:N-1)'*phi1);
% 
% phi2 = 2/N;
% w_test2 = exp(1i*pi*(0:N-1)'*phi2);
% plot(f_fine,abs(w_test1'*A_fine),'LineWidth',1.5)
% hold on
% plot(f_fine,abs(w_test2'*A_fine),'LineWidth',1.5)
%%
% w_combine = exp(-1i*15*pi*1/(N))*0.5*(w_test1 + exp(-1i*2*pi*(N-1)/(N))*w_test2);
% plot(f_fine,abs(w_combine'*A_fine),'LineWidth',1.5)
% w_combine = 0.5*(w_test1 + w_test2);
% plot(f_fine,abs(w_combine'*A_fine),'LineWidth',1.5)
% plot(f_fine,real(w_combine'*A_fine),'LineWidth',1.5)
% ylim([-10,10])
% grid on
% hold on
% plot(f_fine,imag(w_combine'*A_fine),'LineWidth',1.5)
% hold off
% 
% %%
% w_combine = 0.5*(w_test1 + w_test2);
% plot3(real(w_combine'*A_fine), imag(w_combine'*A_fine),f_fine,'LineWidth',1.5)
% hold on 
% w_combine = 0.5*(w_test1 + exp(-1i*4*pi*1/(N))*w_test2);
% plot3(real(w_combine'*A_fine), imag(w_combine'*A_fine),f_fine,'LineWidth',1.5)
% view(2)
% grid on
% hold off


%%
% phi1 = 0;
% w_test1 =  exp(-1i*pi/2)*exp(1i*pi*(0:N-1)'*phi1);
% w_combine = w_test1;
% plot(f_fine,real(w_combine'*A_fine),'LineWidth',1.5)
% hold on 
% plot(f_fine,imag(w_combine'*A_fine),'LineWidth',1.5)
% grid on
% hold off



