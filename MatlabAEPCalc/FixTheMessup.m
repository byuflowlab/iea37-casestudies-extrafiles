% Nicholas F. Baker
% To fix the messup of printing turbine locations in the wrong order
clear, close all;

% Read .csv file
fn = 'FarmList3rdTest9.csv';
LocsIn = csvread(fn);                       % Read in our csv
LocsArray = LocsIn(:);                      % Do so in descending order, down 1st column
LocsArray = LocsArray(1:18000);             % Chop trailing zeros put in from csvread
%LocsArray = LocsArray';                     % Make it a single row
LocsFixed = reshape(LocsArray,[18,1000]);   % Chop the columns
LocsFixed = LocsFixed';                     % Put in the right order
csvwrite('tester.csv',LocsFixed)

% Write to other file
% Row over