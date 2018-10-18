%% Question 1
T = -100:99; % time axis from -100 to 99 ms in steps of 1 ms
lambda = 0.8; % keep lambda the same
N = 20; % keep the number of trials the same
Rand = rand(length(T), N); % Rand mx that is 200 x 20 filled w/ random #s in [0,1]
Spikes = zeros(length(T), N); % Spikes mx that is 200 x 20 initialized with 0s
Spikes(Rand > lambda) = 1; % For all entries of Rand greater than lambda, set the
                           % corresponding indices in Spikes equal to 1

Psth = sum(Spikes'); % Sum up the 1s along all trials for each 1 ms increment
figure; subplot(2,1,1); bar(T,Psth); % In upper subplot, make a bar graph of time vs. Psth;
xlabel('Time in milliseconds')
ylabel('Total number of spikes in N=20 trials')
title('PSTH')
axis tight

subplot(2,1,2); % Lower subplot
[a,b] = find(Spikes ==1); % Find all entries in Spikes equal to 1;
                          % store row indices in 'a' and column in 'b'
a = a - 101; % Shift all elements of 'a' (a vector of row indices)
             % to match the time axis from -100 to 99 ms
scatter(a,b,'.') % Plot - draw dots wherever the coordinate matches (a(i),b(i))
xlabel('Time in milliseconds')
ylabel('Trial number')
title('Raster Plot')

%% Question 2
hold off
clear all

T = 200;
N = 50; % change N=20 to N=50 trials
lambda = 0.1:0.1:0.9; % define a lambda vector to use as an axis in the final plot
Rand = rand(T,N); % Rand mx T x N of random numbers

i = 1; % initialize variable i to 1 for counting
std_axis = zeros(length(lambda),1); % define a s.d. axis of the same dim. as lambda
for i=1:length(lambda) % traverse every element of lambda
  Spikes = zeros(T,N); % define Spikes mx of same dimensions as Rand
  Spikes(Rand>lambda(i)) = 1; % for all entries of Rand above threshold, set
                              % corresponding indices of Spikes equal to 1
  tot_spikes = sum(Spikes); % tot_spikes stores the total # spikes for trials N = 1 to 50
  std_for_this_lambda = std(tot_spikes); % standard deviation in total # spikes
                                         % across trials for current lambda
  std_axis(i) = std_for_this_lambda; % store the s.d. obtained in ith element of std_axis
end

figure
plot(lambda',std_axis) % Plot lambda vector against standard deviation vector
xlabel('Lambda (threshold value)')
ylabel('Standard deviation in total number of spikes across N=50 trials')
title('Standard deviation in total number of spikes per trial as a function of the threshold value')
axis tight

%% Question 3
clear all
hold off

N = 50; % Keep the same as Q2
T = -100:99;

lambda = zeros(length(T), 1); % Define lambda as a companion vector to the
                              % time axis and initialize w/ 0s
for ii = 1:110 % For -100 to +9 ms, set the corresponding entries in lambda to 0.8
  lambda(ii) = 0.8;
end
for jj = 111:170 % For +10 to +70 ms, set lambda to 0.5
  lambda(jj) = 0.5;
end
for kk = 171:200 % For +71 to +200 ms, set lambda to 0.8
  lambda(kk) = 0.8;
end

Rand = rand(length(T),N); % Rand mx 200 x 50 of random numbers

figure
hold on

% For the following lines, each iteration of the for loop corresponds to the current value of lambda
% for the given instance of time. And each iteration creates a scatter plot across the N=1 to 50 trials
% for that time. The 200 vertical scatter plots are consolidated into one to show how the raster plot
% changes with a momentary decrease in lambda from +10 to +70 ms.
for i = 1:1:200
  Spikes = zeros(length(T), N);
  Spikes(Rand>lambda(i)) = 1;
  b = find(Spikes(i,:) == 1); % For the current row, find the indices of all elements
                             % equal to 1 and store in 'b'
  a = zeros(1,length(b)); % Initialize 'a' (with 0s) to be the same dimensions as 'b'
  for iii = 1:length(a) % Traverse every entry of 'a'
    a(iii) = i - 101; % Set all entries to the current row index minus 101
                      % (to account for the shifted time axis, -100 to 99 ms)
  end
  scatter(a,b,'.') % Plot the scatter plot for the given coordinate of time
end

xlabel('Time in milliseconds')
ylabel('Trial number')
title('Raster plot simulation of a neuron response to a stimulus')
axis tight
