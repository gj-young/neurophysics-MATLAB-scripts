clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Question 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

CellCount = 64;
Discretization = 10000;
Angles = linspace(-pi,pi,Discretization);
dt = 1/1000;
wave = 1 + 9 * sin(linspace(0,pi,1500));

PresentationDuration = 5:5:75;     % For incrementing the presentation duration
All_Indices(length(PresentationDuration),1) = 0; % Initalize a vector to store
% the index corresponding to max of Decoding for each presentation duration

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Encoding and decoding for varying %%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% presentation duration times %%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

TuningCurves(1:CellCount, 1:10000) = 1;
CenterPoints = round(linspace(1000,9000,CellCount));

for n=1:CellCount
  TuningCurves(n, CenterPoints(n)-750:CenterPoints(n)+749) = wave;
end

N = 1;     % Set the number of trials equal to 1

for i=1:length(PresentationDuration)
  clear n RandVector Spikes StimInd PopulationVector StimulusSpace m
  clear TestingVectors Decoding
  % Clear variables so don't have to make new ones
  % Only need value of Index for each sampling interval, which will be stored
  % in All_indices vector

  for n = 1:CellCount
    RandVector = rand(Discretization,N);
    Spikes(:,:,n) = (RandVector < (repmat(TuningCurves(n,:)',1,N)) * dt);
  end

  StimInd = 5000;
  for n = 1:CellCount
    PopulationVector(n) = sum(sum(Spikes((StimInd-PresentationDuration(i)):(StimInd+PresentationDuration(i)),:,n)));
  end

  StimulusSpace = 100:100:10000;

  for m=1:length(StimulusSpace)
    for n=1:CellCount
      TestingVectors(n,m) = TuningCurves(n, StimulusSpace(m));
    end
    TestingVectors(:,m) = TestingVectors(:,m)/sum(TestingVectors(:,m));
  end

  for m=1:length(StimulusSpace)
    Decoding(m)=dot(PopulationVector,TestingVectors(:,m));
  end

  [value,Index] = max(Decoding(1,:));
  All_Indices(i) = Index;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Decoding error calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Error(length(PresentationDuration),1) = 0;
% Initialize the vector Error, which will store the errors for using varying
% presentation times (or sampling intervals) in the decoding algorithm

for i=1:length(Error)
  Error(i) = abs(Angles(StimInd) - Angles(StimulusSpace(All_Indices(i))));
end

figure
plot(PresentationDuration,Error, '-ok','MarkerSize',3);
xlabel('Presentation duration in ms')
ylabel('Error in radians')
xlim([5 75])
title({'Error in non-parametric decoding of the population vector',
    'as a function of the sampling interval'})
