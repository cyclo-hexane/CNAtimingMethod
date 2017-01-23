% analyse each set's CNV data

bootstrapIterations = 100;

%%%%%%%%%%% UPDATE PATH TO LOCAL ENVIRONMENT %%%%%%%%%%%%

root = '/Users/danieltemko/Google Drive/evocolab/scripts/cnv_timing/model_files/data/';
setNames = {'Polyp.08.WGS', 'Set.02', 'Set.03', 'Set.04', 'Set.05', 'Set.09.Distal', 'Set.09.Proximal', 'Set.10'};

for j = 1:length(setNames)
    set = setNames{j}
    
    %make the directory for results outputs
    dir = strcat(root, 'model_results/', set);
    mkdir(dir);
    
    % read in mutation information after processing in R
    file = strcat(root, 'matlab/', set, '/', set, '.matlab.input.txt');
    formatSpec = '%s%f%f%f%f%f%f%f';
    data = readtable(file,'Delimiter',',', 'Format', formatSpec);
    points = length(data.alpha);
    % read in total length of the diploid regions and number of mutations
    % in these regions from R
    file = strcat(root, 'matlab/', set, '/', set, '.normal.stats.txt');
    formatSpec = '%f%f';
    diploidStats = readtable(file,'Delimiter',',', 'Format', formatSpec);
    

    %add case numbers to the data, depending on the number of copies of the
    %A allele and B allele
    data.caseNum = repmat(0,points,1);
    for i = 1:points
            data.caseNum(i) = InferCase(data.a2(i), data.b2(i));
    end
    
    %estimate timings and add to the table
    [ t1EstSingle, t1VarSingle, t2EstSingle, t2VarSingle ] = EstimateTimingsSingly( data );
    data.t1_est_single = t1EstSingle';
    data.t2_est_single = t2EstSingle';
    [TEst, t1EstJoint] = EstimateTimings(data, diploidStats.length, diploidStats.mutations);
    data.t1_est_joint = t1EstJoint'; 
    data.T_est_joint = repmat(double(TEst),points,1);
    
    %bootstrapping
    [ t1EstArray  ] = SingleBootstrap( data, diploidStats.length, bootstrapIterations );
    clear t1Diff
    for i = 1:bootstrapIterations
        t1Diff(:,i) = t1EstArray(:,i) - data.t1_est_single;
    end
    t1SquareErrors = t1Diff .* t1Diff;
    clear MSE
    for i = 1:points
        MSE(i) = mean(t1SquareErrors(i,:));
    end
    data.t1_MSE_single = MSE';
    [ TEstArray, t1EstArray  ] = JointBootstrap( data, 't1_est_joint', 'T_est_joint', diploidStats.length, bootstrapIterations );
    clear t1Diff
    for i = 1:bootstrapIterations
        t1Diff(:,i) = t1EstArray(:,i) - data.t1_est_joint;
    end
    
    % plot to check for bias
    h = figure;
    boxplot(t1Diff')
    title(set)
    ylabel('t1 Error')
    file = strcat(root, 'model_results/', set, '/', set, '.t1_errors.pdf');
    print(h,file,'-dpdf')
    
    t1SquareErrors = t1Diff .* t1Diff;
    clear MSE
    for i = 1:points
        MSE(i) = mean(t1SquareErrors(i,:));
    end
    data.t1_MSE_joint = MSE'
    TDiff = TEstArray - data.T_est_joint(1);
    
    % plot to check for bias
    h = figure;
    boxplot(double(TDiff))
    title(set)
    ylabel('T Error')
    file = strcat(root, 'model_results/', set, '/', set, '.T_errors.pdf');
    print(h,file,'-dpdf')
    
    TSquareErrors = TDiff .* TDiff;
    TMSE = mean(TSquareErrors);
    data.T_MSE_joint = repmat(double(TMSE),points,1)
    
    
    % visualise relationship between t1/(t1+t2) single and joint estimates
    data.ratio = data.t1_est_single ./ (data.t1_est_single + data.t2_est_single);
    
    h = figure;
    maxY = max(max(data.t1_est_joint, data.t1_est_single));
    scatter(data.ratio, data.t1_est_joint)
    xlim([0, 1])
    ylim([0, maxY * 1.05]);
    title(set)
    xlabel('ratio')
    ylabel('Joint t1 estimate')
    file = strcat(root, 'model_results/', set, '/', set, '.early.late.ratio.and.joint.estimates.pdf');
    print(h,file,'-dpdf');
    
    % same relationship for single estimates
    h = figure;
    maxY = max(max(data.t1_est_joint, data.t1_est_single));
    scatter(data.ratio, data.t1_est_single)
    xlim([0, 1])
    ylim([0, maxY * 1.05]);
    title(set)
    xlabel('ratio')
    ylabel('Single t1 estimate')
    file = strcat(root, 'model_results/', set, '/', set, '.early.late.ratio.and.single.estimates.pdf');
    print(h,file,'-dpdf');
    
    %visualise the changes in MSE due to joint estimates
    maxEst = max(max(([data.t1_MSE_single, data.t1_MSE_joint])));
    limits = [0 maxEst*1.1];
    h = figure;
    scatter(data.t1_MSE_single, data.t1_MSE_joint)
    xlim(limits);
    ylim(limits);
    title(set)
    xlabel('Single t1 MSE')
    ylabel('Joint t1 MSE')
    file = strcat(root, 'model_results/', set, '/', set, '.MSE.pdf');
    print(h,file,'-dpdf');
    
    % visualise the changes in estimates due to joint clustering
    maxEst = max(max(([data.t1_est_single, t1EstJoint'])));
    limits = [0 maxEst*1.1];
    h = figure;
    scatter(data.t1_est_single, t1EstJoint);
    xlim(limits);
    ylim(limits);
    title(set)
    xlabel('Single t1 Estimates')
    ylabel('Joint t1 Estimates')
    file = strcat(root, 'model_results/', set, '/', set, '.joint.estimates.pdf');
    print(h,file,'-dpdf');
    
    %write results for R analysis
    output = data;
    %delete this column from the output back to R
    output.t2_est_single = [];
    file = strcat(root, 'model_results/', set, '/', set, '.matlab.output.txt');
    writetable(output,file,'Delimiter',' ');
end