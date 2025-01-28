function convergence_test
    % Parameters, realization, n values, epsilon (width of convergence), and K
    M = 500;  nmax = 2000; epsilon = 0.05;  K = 0.5;
    
    fig = figure('Position',[100 100 1000 500], 'Name', 'Convergence Test');
    
    % GUI control config
    uicontrol('Style', 'popupmenu', 'String', {'Normal(2,2)', 'Uniform(2,4)', 'Exponential(λ=2)'}, 'Position', [0 360 120 25], 'Callback', @updateDistribution);
    
    % call function to show different distribution 
    showDistribution('Normal');
    
    function showDistribution(distribution)
        % define RV for different types of distribution
        switch distribution
            case 'Normal'
                actual_mean = 2;
                Xi = normrnd(2, sqrt(2), M, nmax);
            case 'Uniform'
                actual_mean = 3;
                Xi = unifrnd(2, 4, M, nmax);
            case 'Exponential'
                actual_mean = 1/2;
                Xi = exprnd(0.5, M, nmax);     
        end
        
        % Sample mean calulation
        Yj = cumsum(Xi, 2) ./ repmat(1:nmax, M, 1);
        
        % left plot for convergence in probability
        subplot(1,2,1); cla;
        plot(1:nmax,Yj(1:10, :)', 'Color', [0.7 0.7 0.7]);
        hold on;
        
        % boundaries defining epsilon
        yline(actual_mean + epsilon, '--b');
        yline(actual_mean - epsilon, '--b');
        
        % verticle line to count the realization out of boundaries
        vline = xline(1, '-k');
        
        % number of points outside of the band
        out_of_border = Yj(:,1);
        out_idx = abs(out_of_border - actual_mean) > epsilon;

        scatter(ones(sum(out_idx),1), out_of_border(out_idx), 20, 'r', 'filled');
        
        xlabel('n'); ylabel('Y_j');
        dynamictitleHandler = title(sprintf('Sample Paths (p_n = %d/%d paths)', sum(out_idx), M));
        grid on;
        
        % right plot for convergence in probability
        subplot(1,2,2); cla;
        
        % pn calculation
        pn_count = zeros(1, nmax);
        pn = zeros(1, nmax);
    
        for n = 1:nmax
            out_idx = abs(Yj(:,n) - actual_mean) >  epsilon;
            pn_count(n) = sum(out_idx);
            pn(n) = pn_count(n) / M;
        end
        
        % an calculation
        Knmax = floor(K * nmax);
        an = zeros(1, Knmax);
        for n = 1:Knmax
            violation = false(M, 1);

            for i = 1:M
              violation(i) = any(abs(Yj(i,n:end) - actual_mean) > epsilon);
            end
            an(n) = sum(violation) / M;
        end
        
        % plot the pn and an criterion value for convergence in probability
        plot(1:Knmax, an, 'r-', 'DisplayName', 'a_n');
        hold on;
        plot(1:nmax, pn, 'b-', 'DisplayName', 'p_n');
        xlabel('n'); ylabel('pn and an');
        title('Criterion value for convergence in probability (pn and an)');
        legend('show');
        grid on;
        %%%for exponential distribution adjust y limit as mean is different
        if strcmp(distribution, 'Exponential')
            ylim([actual_mean-1, actual_mean+1]);  
        else
            ylim([0 1]);
        end
   
        % Create slider
        delete(findobj(fig, 'Style', 'slider'));
        uicontrol('Style', 'slider', 'Position', [400 20 400 20], 'Min', 1, 'Max', nmax, 'Value', 1, 'SliderStep', [1/nmax, 100/nmax], ...
            'Callback', @(src,~)updateBar(src, Yj, actual_mean, vline, dynamictitleHandler));
    end
    

    %helper functions

    function updateDistribution(src, ~)
        value = src.Value;
        strings = src.String;
        showDistribution(strtok(strings{value}, '('));
    end

    function updateBar(src, Yj, actual_mean, vline, dynamictitleHandler)
        n = round(src.Value);
        vline.Value = n;
        
        % Update scatter points
        subplot(1,2,1);
        delete(findobj(gca, 'Type', 'scatter'));
        out_of_border = Yj(:,n);
        out_idx = abs(out_of_border - actual_mean) > epsilon;
        scatter(n*ones(sum(out_idx),1), out_of_border(out_idx), 20, 'r', 'filled');
        
        % Update title
        set(dynamictitleHandler, 'String', sprintf('Sample Paths (p_n = %d/%d paths outside ε-band)', sum(out_idx), M));
    end
end