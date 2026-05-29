function results = run_logistic(num_classes, num_elements, all_params, classes, class, param_names, glucose, sim_type)
    % Iterate through each class
    for c = 1:num_classes
        current_class = classes{c};
        predictors = []; % To store combined predictor values
        responses = []; % To store combined response values
        
        % Combine data
        for i = 1:num_elements
            % Get predictor values (assuming predictors are in 'low' substructure)
            predictor = all_params{i}';
            
            % Get corresponding class values
            response = strcmp(class{i}, current_class); % Binary (1 if matches current class)
            
            % Append predictors and responses
            predictors = [predictors; predictor];
            responses = [responses; response];
        end
        
        % Create a table for logistic regression
        tbl = array2table(predictors, 'VariableNames', param_names); % Add predictor variables
        tbl.Responses = responses; % Add response variable
        invariant_predictors = {}; % To store names of invariant predictors
        for i = 1:numel(param_names)
            predictor_values = tbl.(param_names{i}); % Extract the column for the current predictor
            % Check if the predictor is invariant
            if numel(unique(predictor_values)) == 1
                invariant_predictors{end+1} = param_names{i}; % Store the name
            end
        end
        
        % Remove invariant predictors from the table
        if ~isempty(invariant_predictors)
            tbl(:, invariant_predictors) = []; % Remove the columns
            table_names = setdiff(param_names, invariant_predictors, 'stable'); % Update param_names
            fprintf('Removed invariant predictors: %s\n', strjoin(invariant_predictors, ', '));
        else
            disp('No invariant predictors found.');
        end

        % Perform logistic regression
        predictor_formula = strjoin(table_names, ' + '); % Join all parameter names with '+'
        formula = sprintf('Responses ~ %s', predictor_formula); % Create the formula string
        glm_model = fitglm(tbl, formula, 'Distribution', 'binomial', 'Link', 'logit');
        
        % Extract significant predictors and p-values
        coefficients = glm_model.Coefficients; % Coefficient table
        significant_idx = find(coefficients.pValue < 0.05 & ~strcmp(coefficients.Row, '(Intercept)')); % Exclude intercept
        significant_predictors = coefficients.Row(significant_idx); % Names of significant predictors
        significant_pvalues = coefficients.pValue(significant_idx); % Corresponding p-values
        significant_odds_ratios = exp(coefficients.Estimate(significant_idx)); % odds ratios for each predictor
        
        % Handle Inf values
        finite_vals = significant_odds_ratios(~isinf(significant_odds_ratios)); % Extract finite values
        if isempty(finite_vals) % Edge case: all are Inf
            significant_odds_ratios(~isfinite(significant_odds_ratios)) = 1; % Replace Inf with 1
        else
            max_finite = max(finite_vals);
            significant_odds_ratios(isinf(significant_odds_ratios)) = max_finite * 1000000; % Replace Inf with slightly higher than max
        end
        
        % Min-max normalization with protection against division by zero
        min_val = min(significant_odds_ratios);
        max_val = max(significant_odds_ratios);
        
        if max_val == min_val
            odds_ratios_norm = ones(size(significant_odds_ratios)); % If all values are the same, set to 1
        else
            odds_ratios_norm = (significant_odds_ratios - min_val) / (max_val - min_val);
        end

        
        % Store the model and significant predictors
        results{c} = struct('Class', current_class, ...
                            'Model', glm_model, ...
                            'SignificantPredictors', significant_predictors, ...
                            'SignificantPValues', significant_pvalues,...
                            'OddsRatios',significant_odds_ratios);
        
        % Display summary for the class
        fprintf('Logistic Regression for Class: %s; %s; %s\n', current_class, sim_type, glucose);
        disp(glm_model);
        
        % Display significant predictors
        if ~isempty(significant_predictors)
            fprintf('Significant Predictors for Class %s:\n', current_class);
            for p = 1:length(significant_predictors)
                fprintf('  %s (p = %.4f); Odds ratio = %.6f \n', significant_predictors{p}, significant_pvalues(p), significant_odds_ratios(p));
            end
        else
            fprintf('No significant predictors for Class %s.\n', current_class);
        end
    end
end