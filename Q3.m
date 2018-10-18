%% Q1
clear all

TimeVector = -10000:1:9999; % Create axis of 20,000 time stamps in 1 ms increments

x = linspace(-pi,pi,20000);
FiringRate = 9 * cos(x) + 1;  % Initialize FiringRate vector as a cosine
FiringRate(1,1:5000) = 1;   % Set the first 5 s to 1 Hz
FiringRate(1,15001:20000) = 1;    % Set the last 5 s to 1 Hz

plot(TimeVector,FiringRate)
xlabel('Time in ms')
ylabel('Firing Rate in Hz')
ylim([0 10])


%% Q2

N = 50;   % Number of trials is 50

dt = 1/1000;

RandMx = rand(length(FiringRate), N);   % 20,000 x 50 initialized with random #s
Spikes = zeros(length(FiringRate), N);  % 20,000 x 50 initialized with 0s
EffectiveFR = FiringRate * dt;   % EffectiveFR between 0 and 1 will be our lambda

for i = 1:length(FiringRate)    % Iterate for each ms increment of time
  temp = RandMx(i,:) < EffectiveFR(i);
% Return all trial numbers (column indices) where the random number is less than
% EffectiveFR for that instant of time and store them in a temporary vector
  Spikes(i,temp) = 1;   % Set the corresponding column indices in Spikes equal to 1
end

figure
PSTH = sum(Spikes');  % Sum the total # spikes for every instant in time
bar(TimeVector,PSTH)  % Create PSTH graph
title('PSTH for sinusoidal time-varying firing rate')
xlabel('Time in ms')
ylabel('Total number of spikes across N=50 trials')

figure
[a,b] = find(Spikes==1);  % Store all indices where there are spikes in [a,b]
a = a - 10000;   % Shift the time axis to desired range -10,000 to 9999 ms
scatter(a,b,'.')
title('Raster plot for time-varying sinusoidal firing rate')
xlabel('Time in ms')
ylabel('Trial number')

%% Q3

 meanISI_First = zeros(1,N);     % Initialize mean, sd, CV vectors
 stdISI_First = zeros(1,N);
 CV_First = zeros(1,N);
 meanISI_Middle = zeros(1,N);
 stdISI_Middle = zeros(1,N);
 CV_Middle = zeros(1,N);

for j = 1:N
  clear ISI_First;
  clear Spikes_current_trial;

  Spikes_current_trial = Spikes(1:5000,j); % Submatrix of spikes: first 5s, jth trial

  ISI_First = diff( find(Spikes_current_trial) ); % Compute all ISIs for first 5 s

  meanISI_First(1,j) = mean(ISI_First); % This vector stores mean ISI for each trial
  stdISI_First(1,j) = std(ISI_First); % This vector stores std in ISI for each trial

  CV_First(j) = stdISI_First(1,j) / meanISI_First(1,j); % Stores CV for each trial
end

for k = 1:N
  clear ISI_Middle;
  clear Spikes_current_trial;

  Spikes_current_trial = Spikes(5001:15000,k); % Submatrix of spikes: middle 10s, kth trial

  ISI_Middle = diff( find(Spikes_current_trial) ); % Compute all ISIs for middle 10 s

  meanISI_Middle(1,k) = mean(ISI_Middle); % This vector stores mean ISI for each trial
  stdISI_Middle(1,k) = std(ISI_Middle); % This vector stores std in ISI for each trial

  CV_Middle(k) = stdISI_Middle(1,k) / meanISI_Middle(1,k); % Stores CV for each trial
end

figure
scatter(meanISI_First,stdISI_First)  % Plot mean vs std in ISIs (first 5 s) for each trial
xlabel('Mean ISI in ms')
ylabel('Standard deviation in ISI in ms')
title('Coefficient of variation values for N=50 trials for first 5 sec.')

figure
scatter(meanISI_Middle,stdISI_Middle)  % Plot mean vs std in ISIs (middle 10 s) for each trial
xlabel('Mean ISI in ms')
ylabel('Standard deviation in ISI in ms')
title('Coefficient of variation values for N=50 trials for middle 10 sec.')


%% Q4

SpikesVector = reshape(Spikes,1000000,1); % Make a spike train for 1000 s
[cor,lags] = xcorr(SpikesVector,'coeff'); % Store time lag in lags and correlation in cor
middleInd = length(SpikesVector);
cor(middleInd) = 0;  
figure
% Plot autocorrelation function for -100 to 100 s around zero lag
plot(lags(middleInd-100000:middleInd+100000),cor(middleInd-100000:middleInd+100000))


xlabel('Time lag in ms')
ylabel('Autocorrelation')
title('Autocorrelation function of 1000 s spike train')
