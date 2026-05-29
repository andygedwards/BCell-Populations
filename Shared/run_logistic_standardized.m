function results = run_logistic_standardized(num_classes, num_elements, all_params, classes, class, param_names, glucose, sim_type)
    % Iterate through each class
    for c = 1:num_classes
        current_class = classes{c};
        predictors = []; % To store combined predictor values
        responses = []; % To store combined response values
        
        % Combine data
        for i = 1:num_elements
            predictor = all_params{i}'; % Get predictor values
            response = strcmp(class{i}, current_class); % Binary (1 if matches current class)
            
            % Append predictors and responses
            predictors = [predictors; predictor];
            responses = [responses; response];
        end
        
        % Standardize predictors (z-score transformation)
        mean_vals = mean(predictors, 1);
        std_vals = std(predictors, 0, 1);
        std_vals(std_vals == 0) = 1; % Avoid division by zero for invariant predictors
        standardized_predictors = (predictors - mean_vals) ./ std_vals;
        
        % Create table for logistic regression
        tbl = array2table(standardized_predictors, 'VariableNames', param_names); 
        tbl.Responses = responses; 
        
        % Identify and remove invariant predictors
        invariant_predictors = {};
        for i = 1:numel(param_names)
            predictor_values = tbl.(param_names{i});
            if numel(unique(predictor_values)) == 1
                invariant_predictors{end+1} = param_names{i}; % Store name
            end
        end
        
        % Remove invariant predictors
        if ~isempty(invariant_predictors)
            tbl(:, invariant_predictors) = []; % Remove columns
            table_names = setdiff(param_names, invariant_predictors, 'stable'); % Update param_names
            fprintf('Removed invariant predictors: %s\n', strjoin(invariant_predictors, ', '));
        else
            disp('No invariant predictors found.');
        end
        
        % % Compute and add quadratic terms
        % quadratic_names = strcat(table_names, '_squared');
        % for i = 1:numel(table_names)
        %     tbl.(quadratic_names{i}) = tbl.(table_names{i}) .^ 2;
        % end

        % Perform logistic regression with quadratic terms
        predictor_formula = strjoin(table_names, ' + '); % strjoin([table_names, quadratic_names], ' + ');
        formula = sprintf('Responses ~ %s', predictor_formula);
        glm_model = fitglm(tbl, formula, 'Distribution', 'binomial', 'Link', 'logit');
        
        % Extract standardized coefficients
        coefficients = glm_model.Coefficients;
        significant_idx = find(coefficients.pValue < 0.05 & ~strcmp(coefficients.Row, '(Intercept)'));
        significant_predictors = coefficients.Row(significant_idx);
        significant_pvalues = coefficients.pValue(significant_idx);
        significant_betas = coefficients.Estimate(significant_idx); % Standardized coefficients
        
        % Compute predicted probabilities
        predicted_probs = glm_model.Fitted.Response; 
        
        % Store results, including predicted probabilities and standardized predictors
        results{c} = struct('Class', current_class, ...
                            'Model', glm_model, ...
                            'SignificantPredictors', significant_predictors, ...
                            'SignificantPValues', significant_pvalues, ...
                            'StandardizedCoefficients', significant_betas, ...
                            'PredictedProbabilities', predicted_probs, ...
                            'StandardizedPredictors', standardized_predictors);
                            % 'QuadraticPredictors', tbl(:, quadratic_names)); % Store quadratic predictor values

        % Display summary
        fprintf('Logistic Regression for Class: %s; %s; %s\n', current_class, sim_type, glucose);
        disp(glm_model);
    end
end
