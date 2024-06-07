function Form = formularset(IEEEstd738Para)
    % Steady-state heat balance 
    Form.I = @(q_c, q_r, q_s, R) sqrt((q_c + q_r - q_s) / R);
    % Non-steady-state heat balance 
    Form.DT_C = @(Dt, I, R, mC_p, q_c, q_r, q_s) (q_s - q_c - q_r + (I)^2 * R) / mC_p * Dt;
    % Convective heat loss
    Form.N_Re = @(D_0, V_w, mu_f, rho_f) D_0 * rho_f * V_w / mu_f;
    % Temperature of the boundary layer
    Form.T_film = @(T_a, T_s) (T_s + T_a) / 2;
    % Forced convection 
    Form.q_c1 = @(K_angle, N_Re, T_a, T_s, k_f) K_angle * (1.01 + 1.35 * N_Re^(0.52)) ...
        * k_f * (T_s - T_a);
    Form.q_c2 = @(K_angle, N_Re, T_a, T_s, k_f) K_angle * (1.01 + 1.35 * N_Re^(0.52)) ...
        * k_f * (T_s - T_a);
    % Wind direction
    Form.K_angle = @(beta) 1.194 - sind(beta) - 0.194 * cosd(2 * beta) + 0.368 * sind(2 * beta);
    Form.q_r = @(D_0, T_a, T_s, epsilon) (17.8 * (D_0) * epsilon ...
        * (((T_s) + 273) / 100)^4 - (((T_a) + 273) / 100)^4);
    % Rate of solar heat gain
    Form.q_s = @(A_prime, Q_se, alpha, theta) alpha * Q_se * sind(theta) * A_prime;
    Form.theta = @(H_c, Z_1, Z_c) acosd(cosd(H_c) * cosd((Z_c - Z_1)));
    % Conductor electrical resistance
    Form.R = @(R_T_low, R_T_high, T_avg, T_high, T_low) ...
        (R_T_high - R_T_low) / (T_high - T_low) * (T_avg - T_low) + R_T_low;
    % Radial temperature gradient within the conductor
    Form.T_s = @(R, T_avg, k_th) T_avg - 1 / 2 * R / (4 * pi * k_th);
    % Conductor heat capacity 
    Form.mC_p = @(m_i, C_p_i) sum(m_i * C_p_i);
    % Dynamic viscosity of air 
    Form.mu_f = @(T_film) 1.458 * 1e-6 * (T_film + 273)^1.5 ...
        / (T_film + 383.4);
    % Air density 
    Form.rho_f = @(H_e, T_film) ...
        (1.293 - 1.525 * 10^-6 * (H_e) + 6.379 * 10^-9 * (H_e)^2) ...
        / (1 + 0.00367 * T_film);
    % Thermal conductivity of air 
    Form.k_f = @(T_film) (2.424 * 1e-2 + 7.477 * 1e-5 * (T_film) ...
        - 4.407 * 1e-9 * (T_film)^2);
    % Altitude of the sun 
    Form.H_c = @(Lat, delta, omega) asind(cosd(Lat) * cosd(delta) * cosd(omega) + sind(Lat) * sind(delta));
    % Solar declination 
    Form.delta = @(N) 23.46 * sind((284 + N) / 365 * 360);
    % Azimuth of the sun 
    Form.Z_c = @(C, chi) C + atand(chi);
    Form.chi = @(Lat, delta, omega) ...
        sind(omega) / (sind(Lat) * cosd(omega) - cosd(Lat) * tand(delta));
    % Total heat flux density at sea level versus 
    Form.Q_s = @(H_c) (IEEEstd738Para.a + IEEEstd738Para.b * (H_c) + IEEEstd738Para.c * (H_c)^2 ...
        + IEEEstd738Para.d * (H_c)^3 + IEEEstd738Para.e * (H_c)^4 ...
        + IEEEstd738Para.f * (H_c)^5 + IEEEstd738Para.g * (H_c)^6);
    % Elevation correction factor 
    Form.Q_se = @(K_solar, Q_s) K_solar * Q_s;
    Form.K_solar = @(H_e) 1 + 1.148 * 1e-4 * H_e + (-1.108 * 1e-8) * (H_e)^2;

end
