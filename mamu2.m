To implement time as a recognizable variable in the program, you can modify the code as follows:


clear; clc

%%
tStep = 1; % Simulation time step [seconds]
% Start and End simulation time
sTime = datetime(2021, 1, 1, 0, 0, 0); 
eTime = datetime(2021, 1, 4, 0, 0, 0);
nTime = seconds(eTime - sTime) / tStep + 1;
lTime = sTime:seconds(tStep):eTime;
rTime = timerange(sTime, eTime, "closed");
%% Input Data
dataTbl = readtable('data.csv');
dataTbl.Time = datetime(dataTbl.YEAR, dataTbl.MO, dataTbl.DY, dataTbl.HR, 0, 0);
dataTT = table2timetable(dataTbl, 'RowTimes', 'Time');
conductorAngle = 50; % 0 degree: direction of North to South
dataTT.Properties.VariableNames = {'T_a', 'V_w', 'W_d', 'solar_radiation'};
% Angle of wind axis calc (-90 to 90 degree)
W_axis = dataTT.W_d;
W_axis(dataTT.W_d > 90) = W_axis(dataTT.W_d > 90) - 180;
W_axis(dataTT.W_d < -90) = W_axis(dataTT.W_d < -90) + 180;
% Beta calc
beta = min(abs((conductorAngle - 90) - W_axis), abs((conductorAngle + 90) - W_axis));
dataTT.beta = beta;
% Initial current
init_I = 1200; % Maximum current is assumed 1200 [A]
dataTT.I = zeros(size(dataTT, 1), 1);
var = {'T_a', 'V_w', 'beta', 'solar_radiation'};
inputTT = retime(dataTT(rTime, var), 'regular', 'linear', 'TimeStep', seconds(1));
time = inputTT.Time;
%% Parameter setting
IEEEstd738Para = struct;

IEEEstd738Para.D_0 = 28.14;
IEEEstd738Para.T_high = 75;
IEEEstd738Para.T_low = 25;
IEEEstd738Para.R_T_low = 7.283 * 10^-5;
IEEEstd738Para.R_T_high = 8.688 * 10^-5;

IEEEstd738Para.Z_1 = 90;
IEEEstd738Para.Lat = 30;

% Radiation factors:
IEEEstd738Para.epsilon = 0.2;
IEEEstd738Para.alpha = 0.2;

% Clear atmosphere
IEEEstd738Para.a = -42.2391; IEEEstd738Para.b = 63.8044; IEEEstd738Para.c = -1.9220; IEEEstd738Para.d = 3.46921 * 10^-2;
IEEEstd738Para.e = -3.61118 * 10^-4; IEEEstd738Para.f = 1.94318 * 10^-6; IEEEstd738Para.g = -4.07608 * 10^-9;

IEEEstd738Para.H_e = 0;
IEEEstd738Para.D_0 = 0.02814;
IEEEstd738Para.D_c = 0.0104;
IEEEstd738Para.mCp = 1310;  % Specific heat capacity
IEEEstd738Para.Dt = tStep;
IEEEstd738Para.K_angle = 1;
IEEEstd738Para.k_th = 1;  % IEEE std and CIGRE 207
%% Calculation - Calculation by MATLAB Function
inputdata = table2array(removevars(timetable2table(inputTT), 'Time')); % 'T_a', 'V_w', 'beta', 'I'

Res.T_avg = zeros(nTime, 1); % Result line temp
Res.T_avg(1, 1) = 40; % Initial line temp
Res.R = zeros(nTime, 1); % Output Resistance
calculated_I = zeros(nTime, 1); % Calculated current will be stored here

% Input data
Tavg = Res.T_avg; % Average temperature
T_abs = inputdata(:, 1); % Absolute temperature
wind_vel = inputdata(:, 2); % Wind velocity
beta = inputdata(:, 3); % Wind direction
solar_radiation = inputdata(:, 4); % Solar radiation

tic
cnt = 0;
for iTime = 1:nTime
    cnt = cnt + 1;
    [Res.T_avg(cnt, 1), calculated_I(cnt), Res.R(cnt, 1)] = IEEEstd738_(Tavg(max(cnt - 1, 1), 1), T_abs(cnt), wind_vel(cnt), beta(cnt), solar_radiation(cnt), init_I, IEEEstd738Para);
end
disp(['Time = ', num2str(toc)]);

Res.Time = time;
Res.current = calculated_I;
Res(:, 1:2)

%% Plot
figure(1);
yyaxis left; plot(Res.Time, Res.T_avg, 'r'); ylabel('Temperature (°C)');
yyaxis right; plot(Res.Time, calculated_I, 'b'); ylabel('Current (A)'); xlabel('Time');
legend('Line Temp', 'Current');
grid on;

