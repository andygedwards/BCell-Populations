function results = IsletCellStatsTwoGroups(y, class)
    % Test for normality (Anderson-Darling test)
    normality_p = adtest(y);
    
    % Test for homogeneity of variances (Bartlett's test)
    homogeneity_p = vartestn(y, class, 'TestType', 'Bartlett', 'Display', 'off');

    % Check the unique groups
    groups = unique(class);
    if numel(groups) ~= 2
        error('This function is designed for exactly two groups.');
    end
    % Split data into the two groups
    group1 = y(find(strcmp(class,groups{1})));
    group2 = y(find(strcmp(class,groups{2})));

    if normality_p > 0.05 && homogeneity_p > 0.05
        % Data meets assumptions for parametric test
        [H,results.pValue, ~, stats] = ttest2(group1, group2);
        results.test = 'Two-sample t-test';
        results.stats = stats;
    else
        % Data does not meet assumptions, use non-parametric test
        [p,H,stats] = ranksum(group1, group2);
        results.pValue = p;
        results.test = 'Mann-Whitney U test';
        results.stats = stats;
    end
    % Store normality and homogeneity test results
    results.normality_p = normality_p;
    results.homogeneity_p = homogeneity_p;
end
