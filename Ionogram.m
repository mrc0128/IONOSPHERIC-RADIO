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
Time = [year, month, day, hour, 0, 0];

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






