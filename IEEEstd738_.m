function [T_avg, I, R] = IEEEstd738_(T_avg, T_a, V_w, beta, t, Rad_solar, IEEEstd738Para)

    persistent T_s T_film rho_f mu_f k_f K_angle chi C Z_c theta q_s DT_C Form

    Form = formularset(IEEEstd738Para);
    
    % Initial average temperature of the boundary layer
    R = Form.R(IEEEstd738Para.R_T_low, IEEEstd738Para.R_T_high, T_avg, IEEEstd738Para.T_high, IEEEstd738Para.T_low); % Resistance at initial temperature
    
    T_s = Form.T_s(R, T_avg, IEEEstd738Para.k_th); % Conductor surface temperature
    
    %% 0 Calculate days(N) and Solar zenith angle
    N = floor(t / (24 * 3600));
    s = rem(t, 24 * 3600);
    omega = (s - 12 * 3600) / 3600 * 15;
    
    %% 1 Calc air condition
    T_film = Form.T_film(T_a, T_s);
    mu_f = Form.mu_f(T_film);
    rho_f = Form.rho_f(IEEEstd738Para.H_e, T_film);
    N_Re = Form.N_Re(IEEEstd738Para.D_0, V_w, mu_f, rho_f);
    k_f = Form.k_f(T_film);
    K_angle = Form.K_angle(beta);
    
    % q_c
    q_c1 = Form.q_c1(K_angle, N_Re, T_a, T_s, k_f);
    q_c2 = Form.q_c2(K_angle, N_Re, T_a, T_s, k_f);
    q_c = max(q_c1, q_c2);
    
    % q_r
    q_r = Form.q_r(IEEEstd738Para.D_0, T_a, T_s, IEEEstd738Para.epsilon);
    
    % q_s
    delta = Form.delta(N);
    chi = Form.chi(IEEEstd738Para.Lat, delta, omega);
    C = 180 * (chi >= 0 & 0 <= omega & omega < 180) ...
        + 180 * (chi < 0) + 180 * (chi < 0 & 0 <= omega & omega < 180);
    Z_c = Form.Z_c(C, chi);
    H_c = Form.H_c(IEEEstd738Para.Lat, delta, omega);
    H_c(H_c < 0) = 0;
    theta = Form.theta(H_c, IEEEstd738Para.Z_1, Z_c);
    Q_s = Form.Q_s(H_c);
    K_solar = Form.K_solar(IEEEstd738Para.H_e);
    Q_se = Form.Q_se(K_solar, Q_s);
    q_s = Form.q_s(IEEEstd738Para.D_0, Q_se, IEEEstd738Para.alpha, theta);
    q_s(q_s < 0) = 0;
    
    R = Form.R(IEEEstd738Para.R_T_low, IEEEstd738Para.R_T_high, T_avg, IEEEstd738Para.T_high, IEEEstd738Para.T_low); % Resistance
    I = Form.I(q_c, q_r, q_s, R);
    DT_C = Form.DT_C(IEEEstd738Para.Dt, I, R, IEEEstd738Para.mCp, q_c, q_r, q_s);
    T_avg = T_avg + DT_C;
    
    T_s = Form.T_s(R, T_avg, IEEEstd738Para.k_th);
end
