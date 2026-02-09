clear; close all; clc;

%% Paramter setup
% Locations
auburn_lat = 32.6;   % Reciever
auburn_lon = -85.5;

chesap_lat = 36.8;  % transmitter 1
chesap_lon = -76.3;

corpus_lat = 27.8;   % transmitter 2
corpus_lon = -97.4;

% Micah Day
year = 2003;
month = 6;
day = 24;
hour = 12;  % UTC
Time = [year, month, day, hour, 0];

% Frequency range
freq_start = 1;      % MHz
freq_end = 10;       % MHz
freq_step = 0.1;     % MHz

fprintf('Receiver: Auburn, AL (%.2f, %.2f)\n', auburn_lat, auburn_lon);
fprintf('Transmitter 1: Chesapeake, VA (%.2f, %.2f)\n', chesap_lat, chesap_lon);
fprintf('Transmitter 1: Chesapeake, VA (%.2f, %.2f)\n', corpus_lat, corpus_lon);
fprintf('Date and Time: %d-%02d-%02d %02d:00 \n', year, month, day, hour);

% calculate range from auburn to chesap

orig_lat = chesap_lat;
orig_lon = chesap_lon;
[range_ches, bearing_ches] = latlon2raz(auburn_lat, auburn_lon, orig_lat, ...
    orig_lon, 'wgs84');

fprintf('\nChesapeake to Auburn:\n');
fprintf('    Range: %.2f km\n', range_ches/1000);
fprintf('    Bearing: %.2f degrees\n', bearing_ches);

% calculate range from auburn to corpus 
orig_lat = corpus_lat;
orig_lon = corpus_lon;
[range_corpus, bearing_corpus] = latlon2raz(auburn_lat, auburn_lon, orig_lat, ...
    orig_lon, 'wgs84');

fprintf('\nCorpus to Auburn:\n');
fprintf('    Range: %.2f km\n', range_corpus/1000);
fprintf('    Bearing: %.2f degrees\n', bearing_corpus);



R12 = 100; % sunspot 
max_range = range_ches/1000; % num range step
num_range = 201; 
range_inc = max_range/(num_range-1);
start_height = 60;
height_inc = 2;
num_heights = 201; % num of height steps
kp = 4; % kp index 
doppler_flag = 1; % doppler parameter 


% generate ionogram for chesapeake transmitter
tic
[iono_pf_grid_ches, iono_pf_grid_5_ches, collision_freq_ches, irreg_ches, iono_te_grid_ches] = gen_iono_grid_2d(chesap_lat, chesap_lon, R12, Time, ...
    bearing_ches, max_range, num_range, range_inc, ...
    start_height, height_inc, num_heights, kp, doppler_flag);
toc

fprintf('Grid Size: %d heights x %d ranges\n', num_heights, num_range);

% ray trace at 5 MHz

fprintf('\nSingle ray trace for 5MHz\n');

% ray param

test_freq = 5.0; %MHz
test_elev = 10; % degrees (first guess)
tol = [1e-7 0.01 25]; % ODE, target, iterations 
nhops = 1; %hops
ir_flag = 1; % irregularities flag 

%trace ray 

tic

[ray_data, ray_path_data, ray_state_vec] = raytrace_2d(auburn_lat, auburn_lon, test_elev, bearing_ches, test_freq, nhops, tol, ir_flag, iono_pf_grid_ches, iono_pf_grid_5_ches, collision_freq_ches, start_height, height_inc, range_inc, irreg_ches);
toc

fprintf('final ground range: %.2f km\n', ray_data.ground_range);




