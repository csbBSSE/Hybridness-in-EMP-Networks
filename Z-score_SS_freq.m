%Comments%
%This code takes network solution file (which is of the Matrix form) as input. It then performs the Z-score standardization of each 
%variable (columns of the solution matrix), before descritizing the (varible) values to get a matrix of 0's and 1's, where each row 
%represents the steady state. The output is the frequency of each unique steady state in the descritized matrix and is written to the excel file.

%clear all;
clear variables;
clear global;
%close all;
%clc;
format short;

wild = readmatrix('wild_run1.xlsx'); 

%load all files in the folder
files = dir('*.dat');
sorted_files = natsortfiles(files);
%T1 = {};
for i=1:length(sorted_files)  %length
    T1 = importdata(sorted_files(i).name);
   
xx = size(T1);
col_length_T1 = xx(2);

T = T1(:, 4:col_length_T1);

yy = size(T);
col_length_T = yy(2);

L = length(T(:, 1)); %row length of RACIPE solution file

score = zeros(L, col_length_T);
%Calculating Z-Score Matrix. Formaula is (datapoint-mean/std)%
for j = 1:col_length_T  % Taking first through last columns %for grhl network's last column is snail that should be dropped so that j = 1:col_length_T-1 
    for k = 1:L           % Taking first through last rows
    score(k, j) = (T(k, j) - mean(T(:, j)))/std(T(:, j)); % dividing a vector by a number (i.e. std)
    end
end

Z = score > zeros(size(score)); %Changing entries into one's and zero's 
% Extracting Unique rows 
[ii,jj,kk] = unique(Z,'rows','stable'); % 'stable' just keeps the order of the rows as is in the original matrix ('Z' in this case)
frequency = histc(kk,1:numel(jj)); % Calculating Frequency
OverallFrequency = [ii frequency];
sort_freq = sortrows(OverallFrequency, 1:col_length_T); %(1:22); % (-4 indicates descending order based on column 4)
%writematrix(sort_freq, 'wild_run1.xlsx')

%replacing 0's by -1 for frust calc

%......A single command that replaces 0's by -1 would be network_SS(network_SS ==0 ) = -1 or 2*network_SS - 1

%network_SS = sortrows(ii, 1:col_length_T);
%network_SS_1 = network_SS-1; 
%converted_states = network_SS_1 + network_SS; 
%length(converted_states);


%wild_type_network data made global variable % see at the top of the
%program
mutant = sort_freq;  % data obtained after running simulation of mutant network

wild_rows_not_in_mutant = find(~ismember(wild(:, 1:57), mutant(:, 1:57), 'rows')); % finding index of those rows in wild that are not in mutant
b = wild(wild_rows_not_in_mutant, 1:57); %actual rows in wild that are not in mutant
c = zeros(1, size(b, 1));   %size(b, 1) counts number of rows in b
d = [b c']; %rows (with zero frequencies) in wild that are not in mutant
ad = [mutant; d]; % all rows (with frequencies) in %mutant 58plus those that are in wild but not in mutant % wild that are in mutant as well as those that are not in mutant
mutant_complete = sortrows(ad, 1:57); %sorting rows based on frequencies
mutant_frequencies = mutant_complete(:, 58); %mut_rows_in_wild_with_freq(:, 5);  % extracting frequency column of mutant

mut_rows_not_in_wild = find(~ismember(mutant(:, 1:57), wild(:, 1:57), 'rows'));
b_prime = mutant(mut_rows_not_in_wild, 1:57);
c_prime = zeros(1, size(b_prime, 1));
d_prime = [b_prime c_prime'];
ad_prime = [wild; d_prime]; % all rows (with frequencies) in mutant that are in mutant as well as those that are not in mutant
wild_complete = sortrows(ad_prime, 1:57);
wild_frequencies = wild_complete(:, 58);

%xlswrite('ZEB-ZEB_1-2.xlsx', wild_complete) %writing altered states of wild-network (with frequencies) to excel file
%xlswrite('ZEB-ZEB_1-2.xlsx', mutant_frequencies, 'Sheet1', 'I') % writing mutant-frequencies in the mut column of the above excel file
%-----Assigning names to 23rd and 24th columns containing wild & mutant frequencies
%col_header(1,21:22) = {'wild', 'mut'}; % Assigning names to 23rd and 24th columns, containing wild & mutant frequencies, as wild and mutant
%xlswrite('ZEB-ZEB_1-2.xlsx',col_header); %naming columns of the excel file by the above created names

filename = sorted_files(i).name;
fullfilename = strcat(filename(1:length(filename)-13),'.xlsx');

%writematrix(converted_states, fullfilename);

writematrix(wild_complete, fullfilename);
writematrix(mutant_frequencies, fullfilename, 'Sheet', 1, 'Range','BG1:BG70000');   % change V to X for EMT_files

%col_header(1,21:22) = {'wild', 'mut'};
%writematrix(col_header, fullfilename);

end
