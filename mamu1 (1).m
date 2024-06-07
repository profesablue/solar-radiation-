To implement and change time to be a recognizable variable in the program, you can modify the code as follows:


clear; clc

%%
tStep = 3600; % Simulation time step [seconds]
% Start and End simulation time
sTime = datetime(2021,1,1,0,0,0); 
eTime = datetime(2021,1,2,0,0,0);
lTime = sTime:tStep:eTime;
rTime = timerange(sTime,eTime,"closed");
%% Input Data
dataTbl = readtable('data.csv');

% Extract time columns and convert to datetime
timeColumns = dataTbl(:, 1:4);
timeData = table2array(timeColumns);
timeArray = datetime(timeData,'InputFormat','yyyy-MM-dd HH:mm:ss');

% Create timetable with time as row times
dataTT = array2timetable(table2array(dataTbl(:, 5:end)), 'RowTimes', timeArray);
conductorAngle = 50; % 0 degree: direction of North to South
dataTT.Properties.VariableNames = {'Solar_Radiation', 'T_a', 'V_w', 'W_d'};
% angle of wind axis calc (-90 to 90 degree)
W_axis = dataTT.W_d;
W_axis(dataTT.W_d > 90) = W_axis(dataTT.W_d > 90) - 180;
W_axis(dataTT.W_d < -90) = W_axis(dataTT.W_d < -90) + 180;
% beta calc
beta = min(abs((conductorAngle - 90) - W_axis), abs((conductorAngle + 90) - W_axis));
dataTT.beta = beta;
% initial current
init_I = zeros(size(dataTT, 1), 1);
dataTT.I = init_I;
var = {'T_a', 'V_w', 'beta', 'Solar_Radiation'};
inputTT = retime(dataTT(rTime, var), 'regular', 'linear', 'TimeStep', seconds(1));
t = seconds(inputTT.Time - inputTT.Time(1)); % Time in seconds

% Define the time variable
time = lTime;

%% Parameter setting
IEEEstd738Para = struct;

IEEEstd738Para.D_0 = 28.14;
IEEEstd738Para.T_high = 75;
IEEEstd738Para.T_low = 25;
IEEEstd738Para.R_T_low = 7.283 * 10^-5;
IEEEstd738Para.R_T_high = 8.688 * 10^-5;

IEEEstd738Para.Z_1 = 90;
IEEEstd738Para.Lat = 30;

% radiation factors:
IEEEstd738Para.epsilon = 0.2;
IEEEstd738Para.alpha = 0.2;

%% Clear atmosphere parameters
IEEEstd738Para.a = -42.2391; IEEEstd738Para.b = 63.8044; IEEEstd738Para.c = -1.9220; IEEEstd738Para.d = 3.46921 * 10^-2;
IEEEstd738Para.e = -3.61118 * 10^-4; IEEEstd738Para.f = 1.94318 * 10^-6; IEEEstd738Para.g = -4.07608 * 10^-9;

IEEEstd738Para.H_e = 0;
IEEEstd738Para.D_0 = 0.02814;
IEEEstd738Para.D_c = 0.0104;
IEEEstd738Para.mCp = 1310; % specific heat capacity
IEEEstd738Para.Dt = tStep;
IEEEstd738Para.K_angle = 1;
IEEEstd738Para.k_th = 1; % IEEE std and CIGRE 207

%% Calculation - MATLAB Function
inputdata = table2array(removevars(timetable2table(inputTT), 'Time')); % 'T_a', 'V_w', 'beta', 'I'

Res.T_avg = zeros(length(inputdata), 1); % Result line temp
Res.T_avg(1) = 40; % Initial line temp
Res.R = zeros(length(inputdata), 1); % Output Resistance
calculated_I = zeros(length(inputdata), 1); % Calculated current will be stored here

% input data
Tavg = Res.T_avg; % Average temperature
T_abs = inputdata(:, 1); % Absolute temperature
wind_vel = inputdata(:, 2); % Wind velocity
beta = inputdata(:, 3); % Wind direction
solar_radiation = inputdata(:, 4); % Solar radiation

tic
for i = 1:length(inputdata)
    [Tavg(i), calculated_I(i), Res.R(i)] = IEEEstd738_(Tavg(max(i-1,1)), T_abs(i), wind_vel(i), beta(i), t(i), IEEEstd738Para, solar_radiation(i));
end
toc

% output data
% update I
inputTT.I = calculated_I;
outputTT = inputTT;
outputTT.Time = lTime(1:length(outputTT.I));
% update beta
beta(1) = 0;
outputTT.beta = beta;
outputTT.Properties.VariableNames{5} = 'I';
outputTT.Properties.VariableNames{4} = 'Solar_Radiation';
outputTT.Properties.VariableNames{3} = 'Wind_Direction';
outputTT.Properties.VariableNames{2} = 'Wind_Speed';
outputTT.Properties.VariableNames{1} = 'Temperature';

writetimetable(outputTT, 'output.csv');

%% Plot result
figure
yyaxis left
plot(lTime, calculated_I, '-b', 'LineWidth', 2);
xlabel('Time');
ylabel('Calculated Current [A]');
yyaxis right
plot(lTime, Tavg, '-r', 'LineWidth', 2);
ylabel('Average Temperature [C]');
grid on
legend('Calculated Current', 'Average Temperature');
title('Line Current and Average Temperature vs Time');

