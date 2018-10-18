clear all

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Question 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

NeuronCountVector = 4:4:64; % For incrementing the neuron count in the popu vector
Discretization = 10000;
Angles = linspace(-pi,pi,Discretization);
dt = 1/1000;
wave = 1 + 9 * sin(linspace(0,pi,1500));

All_Indices( length(NeuronCountVector), 1) = 0; % Initalize a vector to store
% the index corresponding to max of Decoding for each population vector length

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% Encoding and decoding loop for varying %%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%% population vector lengths %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

for i=1:length(NeuronCountVector)
  clearvars -except Decoding NeuronCountVector Discretization Angles dt wave i All_Indices
% Use the variables defined in the loop through NeuronCountVector to temporarily
% store values for the purpose of finding Index and populating All_Indices
% Clear them after each iteration

  TuningCurves(1:NeuronCountVector(i), 1:10000) = 1;
  CenterPoints = round(linspace(1000,9000,NeuronCountVector(i)));

  for n=1:NeuronCountVector(i)
    TuningCurves(n, CenterPoints(n)-750:CenterPoints(n)+749) = wave;
    end

  N = 1; % Set number of trials equal to 1

  for n = 1:NeuronCountVector(i)
    RandVector = rand(Discretization,N);
    Spikes(:,:,n) = (RandVector < (repmat(TuningCurves(n,:)',1,N)) * dt);
  end

    StimInd = 5000;
    for n = 1:NeuronCountVector(i)
      PopulationVector(n) = sum(sum(Spikes((StimInd-50):(StimInd+50),:,n)));
    end


  StimulusSpace = 100:100:10000;

  for m=1:length(StimulusSpace)
    for n=1:NeuronCountVector(i)
        TestingVectors(n,m) = TuningCurves(n, StimulusSpace(m));
    end
    TestingVectors(:,m) = TestingVectors(:,m)/sum(TestingVectors(:,m));
  end

  for m=1:length(StimulusSpace)
    Decoding(i,m)=dot(PopulationVector,TestingVectors(:,m));
  end

  [value,Index] = max(Decoding(i,:));
  All_Indices(i) = Index;
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%% Decoding error calculation %%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Error(length(NeuronCountVector),1) = 0;
% Initialize the vector Error, which will store the errors for
% using varying lengths of the population vector in the decoding algorithm

for i=1:length(Error)
  Error(i) = abs(Angles(StimInd) - Angles(StimulusSpace(All_Indices(i))));
end

figure
plot(NeuronCountVector,Error, '-ok','MarkerSize',3);
xlabel('Number of cells')
ylabel('Error (radians)')
xlim([4 64])
title({'Error in non-parametric decoding of the population vector as a function',
    'of the number of cells included in the population vector'})
