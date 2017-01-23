%script takes reasonable assumptions for ground truth values of the ti ant
%T, simulates data based on these values and uses the inference pipeline to
%estimates the ti and T. These estimates are then assessed against the
%known ground truth. 
clear all

T = 10; %years
mu = 5e-8; %mutaions per base pair per year

%case 1 only test
clear dummy1
times1 = [5];
Time = 10;
dummy1.CNV_ID = [1];
dummy1.length = [9e7];
dummy1.a2 = [0];
dummy1.b2 = [2];
dummy1.t1 = times1*mu;
points = length(dummy1.a2); 
dummy1.T = repmat(Time*mu ,1, points);
diploidLength = 9e9;

%add case numbers
for i = 1:points
    dummy1.caseNum(i) = InferCase(dummy1.a2(i), dummy1.b2(i));
end

%Generate mutations 
% inputs to this function ar (data, t1ColName, TColName, diploidLength, iterations)
[alpha, beta, diploidRegionMutations] = GenerateMutations(dummy1, 't1', 'T', diploidLength, 1)
dummy1.alpha = alpha;
dummy1.beta = beta;

%estimate the ti and T
[ TEst, t1Est ] = EstimateTimings( dummy1, diploidLength, diploidRegionMutations )

%bootstrapping
iterations = 200;
[ TEstArray, t1EstArray  ] = JointBootstrap( dummy1, 't1', 'T', diploidLength, iterations );
for i = 1:iterations
        t1DummyDiff(:,i) = t1EstArray(:,i) - dummy1.t1;
end
boxplot(t1DummyDiff)    

%case 2 only test
dummy1.a2 = 2;

%case 3 only test
dummy1.a2 = 1;

%add case numbers
for i = 1:points
    dummy1.caseNum(i) = InferCase(dummy3.a2(i), dummy3.b2(i));
end

%Generate mutations and estimate the ti and T 10 times
[ TEstArray, t1EstArray  ] = JointBootstrap(dummy1, 't1', 'T', diploidLength, 10);
mean(TEstArray)
mean(t1EstArray)

%mixed case example
clear data
T = 10;
mu = 5e-8;
times1 = [5 4 7 3];
Time = 10;
dummymixed.CNV_ID = [1 2 3 4];
dummymixed.length = [9e7 2.4e8 1.9e8 2e8];
dummymixed.a2 = [0 1 1 2];
dummymixed.b2 = [3 2 3 2];
dummymixed.t1 = times1*mu;
points = length(dummy1.a2); 
dummymixed.T = repmat(Time*mu ,1, points);
normalLength = 1.5e9;

%print assumption sets to screen
dummy1
dummy2
dummy3
dummymixed

dummy1.t1 = times1*mu;
points = length(dummy1.a2); 
dummy1.T = repmat(Time*mu ,1, points);
normalLength = 1.5e9;

%case 2 only test
clear dummy2
dummy2 = dummy1;
dummy2.a2 = 2;

%case 3 only test
clear dummy3
dummy3 = dummy1;
dummy3.a2 = 1;

%mixed case example
clear data
T = 10;
mu = 5e-8;
times1 = [5 4 7 3];
Time = 10;
dummymixed.CNV_ID = [1 2 3 4];
dummymixed.length = [9e7 2.4e8 1.9e8 2e8];
dummymixed.a2 = [0 1 1 2];
dummymixed.b2 = [3 2 3 2];
dummymixed.t1 = times1*mu;
points = length(dummymixed.a2); 
dummymixed.T = repmat(Time*mu ,1, points);
normalLength = 1.5e9;

%print assumption sets to screen
dummy1
dummy2
dummy3
dummymixed

%test against 

%add case numbers
for i = 1:points
    dummy1.caseNum(i) = InferCase(dummy1.a2(i), dummy1.b2(i));
end

%generate data table
clear data
T = 10;
mu = 5e-8;
times1 = [5];
Time = 10;
dummy1.CNV_ID = [1];
dummy1.length = [9e7];
dummy1.a2 = [1];
dummy1.b2 = [3];
dummy1.t1 = times1*mu;
points = length(dummy1.a2); 
dummy1.T = repmat(Time*mu ,1, points);
normalLength = 0;

[alpha beta normalMutations] = GenerateMutations(dummy1, normalLength, 1);
dummy1.alpha = alpha';
dummy1.beta = beta';

%add case numbers
for i = 1:points
    dummy1.caseNum(i) = InferCase(dummy1.a2(i), dummy1.b2(i));
end

[ TEst, t1Est ] = EstimateTimings( dummy1, normalLength, normalMutations )