function [ t1EstArray  ] = SingleBootstrap( data, normalLength, iterations )

[alphaArray, betaArray, normalMutationsArray] = GenerateMutations(data, 't1_est_single', 'T_est', normalLength, iterations);

for i = 1:iterations
    dummyData = data;
    dummyData.alpha = alphaArray(:,i);
    dummyData.beta = betaArray(:,i);
    dummyNormalMutations = normalMutationsArray(1,i);
    fprintf('Bootstrap Iteration %d...\n',i)
    [ t1EstSingle, t1VarSingle, t2EstSingle, t2VarSingle ] = EstimateTimingsSingly( dummyData );
    t1EstArray(:,i) = t1EstSingle';
end

end

