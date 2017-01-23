function [ t1Est, t1Var, t2Est, t2Var ] = EstimateTimingsSingly( data )

points = length(data.alpha);

%loop through all data rows and extract esimate
for i = 1:points
    copiedAlleles = (data.a2(i) > 1) + (data.b2(i) > 1);
    copiedLength = copiedAlleles * data.length(i);
    lambda1Est = data.alpha(i);
    lambda1Var = data.alpha(i);
    t1Est(i) = lambda1Est / copiedLength;
    t1Var(i) = lambda1Var / copiedLength^2;
    earlyLowVAFEst = ((data.a2(i) == 1) + (data.b2(i) == 1)) * t1Est(i) * data.length(i); 
    lateMutationsEst = data.beta(i) - earlyLowVAFEst;
    remainingLength = (data.a2(i) + data.b2(i)) * data.length(i);
    lambda2Est = lateMutationsEst;
    lambda2Var = lateMutationsEst;
    t2Est(i) = lambda2Est / remainingLength;
    t2Var(i) = lambda2Var / remainingLength^2;
end
end
	


