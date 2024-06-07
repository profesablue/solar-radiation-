To implement or change the time to be a recognizable variable in the program, you can modify the code as follows:


clear; clc

%%
tStep=1; %Simulation time step [seconds]
% Start and End simulation time
sTime=datetime(2021,1,1,3,0,0); 
eTime=datetime(2021,1,4,2,0,0);
nTime=seconds(eTime-sTime)/tStep+1;
lTime=sTime:seconds(tStep):eTime;
rTime=timerange(sTime,eTime,"closed");

% Define the time variable
time = lTime;

%% Input Data
dataTbl = readtable('exampledata.csv');
dataTbl.Time = datetime(dataTbl.Time, 'InputFormat', 'yyyy/MM/dd HH:mm');
dataTT = table2timetable(dataTbl, 'RowTimes', 'Time');
conductorAngle = 50; % 0 degree: direction of North to South
dataTT.Properties.VariableNames = {'T_a','V_w','W_d','SR'};
% angle of wind axis calculation (-90 to 90 degree)
W_axis = dataTT.W_d;
W_axis(dataTT.W_d > 90) = W_axis(dataTT.W_d > 90) - 180;
W_axis(dataTT.W_d < -90) = W_axis(dataTT.W_d < -90) + 180;
% beta calculation
beta = min(abs((conductorAngle - 90) - W_axis), abs((conductorAngle + 90) - W_axis));
dataTT.beta = beta;
% initial current
init_I = zeros(size(dataTT, 1), 1); % Initial current set to zero for solar radiation calculation
dataTT.I = init_I;
var = {'T_a','V_w','beta','SR'};
inputTT = retime(dataTT(rTime, var), 'regular', 'linear', 'TimeStep', seconds(1));
t = second(inputTT.Time, 'secondofday') + 24 * 3600 * day(inputTT.Time, 'dayofyear');
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

%% Clear atmosphere
IEEEstd738Para.a = -42.2391; IEEEstd738Para.b = 63.8044; IEEEstd738Para.c = -1.9220; IEEEstd738Para.d = 3.46921 * 10^-2;
IEEEstd738Para.e = -3.61118 * 10^-4; IEEEstd738Para.f = 1.94318 * 10^-6; IEEEstd738Para.g = -4.07608 * 10^-9;

IEEEstd738Para.H_e = 0 ;
IEEEstd738Para.D_0 = 0.02814;
IEEEstd738Para.D_c = 0.0104;
IEEEstd738Para.mCp = 1310;  % Specific heat capacity
IEEEstd738Para.Dt = tStep;
IEEEstd738Para.K_angle = 1;
IEEEstd738Para.k_th = 1;    % IEEE std and CIGRE 207
%% Calculation - Calculation by MATLAB Function
inputdata = table2array(removevars(timetable2table(inputTT), 'Time')); % 'T_a','V_w','beta','I'

Res.T_avg = zeros(nTime, 1);           % Result line temp
Res.T_avg(1, 1) = 40;                  % Initial line temp
Res.R = zeros(nTime, 1);               % Output Resistance
calculated_I = zeros(nTime, 1);        % Calculated current will be stored here

% Input data
Tavg = Res.T_avg;                      % Average temperature
T_abs = inputdata(:, 1);               % Absolute temperature
wind_vel = inputdata(:, 2);            % Wind velocity
beta = inputdata(:, 3);                % Wind direction
solar_radiation = inputdata(:, 4);     % Solar radiation

tic
cnt = 0;
for iTime = sTime:seconds(tStep):eTime
    cnt = cnt + 1;
    [Res.T_avg(cnt, 1), calculated_I(cnt), Res.R(cnt, 1)] = IEEEstd738_(Tavg(max(cnt-1, 1), 1), T_abs(cnt), wind_vel(cnt), beta(cnt), t(cnt), ...
            IEEEstd738Para, solar_radiation(cnt));
end
toc
%%
current_abs_val = abs(calculated_I);

table_current = array2table([string(time'), current_abs_val], 'VariableNames', {'Time', 'Calculated Current'});
writetable(table_current, 'savedData2.csv');
%% Simulation result

figure=figure();

title('Line temp and conditions over two days in summer')
ax = gca;

yyaxis left
plot(time, current_abs_val, '-r', 'LineWidth', 1.5);
hold on 

ylabel('Current-amps & Solar Heat Intensity-w/m2');

hold off

yyaxis right
[ax.YAxis.Color] = deal([0,0,0],[0,0,0]);

plot(dataTT.Time, dataTT.T_a, '--b', 'LineWidth', 1.5);
hold on 
plot(dataTT.Time, dataTT.V_w, '-g', 'LineWidth', 1.5);

if exist('Res', 'var')
    if length(time) == length(Res.T_avg)
        plot(time', Res.T_avg, '-k', 'LineWidth', 1.5);
    else
        disp('No equal result var length');
    end
    
    annotation('textbox', ...
    [0.34809523809523 0.275476190476192 0.203285714285714 0.200634920634921], ...
    'Color', [0 0 0], ...
    'String', 'Line temp-C (MATLAB)', ...
    'LineStyle', 'none', ...
    'FitBoxToText', 'off');
end

ylabel('Wind Speed-mps & Air Temperature-C');

hold off

annotation('textbox', ...
    [0.429571428571428 0.606349206349207 0.163285714285714 0.120634920634921], ...
    'Color', [1 0 0], ...
    'String', {'Line current'}, ...
    'LineStyle', 'none', ...
    'FitBoxToText', 'off');

annotation('textbox', ...
    [0.429393675741416 0.203581603581606 0.163285714285714 0.120634920634922], ...
    'Color', [0 0 1], ...
    'String', 'Air Temp-C', ...
    'LineStyle', 'none', ...
    'FitBoxToText', 'off');

annotation('textbox', ...
    [0.437907100757787 0.0611314611314635 0.163285714285713 0.120634920634922], ...
    'Color', [0.3 1.0 0.3], ...
    'String', 'Wind spd', ...
    'LineStyle', 'none', ...
    'FitBoxToText', 'off');

ax.XLim = [time(1) time(end)];



