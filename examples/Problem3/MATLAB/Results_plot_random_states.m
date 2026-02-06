clear;
clc;

% --- First set (N=3) ---
std_d = [];
mean_d = [];
sigma_d = [];
mean_delta_d = [];
dims = 4:4;
mean_ab_d=[];
for d = dims
    filename = sprintf('random_state_d%d_N4_r%d.mat', d, d-1);
    load(filename, 'test_results');    % Make sure test_results is in the .mat file
    N = size(test_results, 1);
    
    % Compute mean
    delta = (test_results(:,2)-test_results(:,1))./(1-test_results(:,1));
    meanVals = mean(delta);
    stdVals = std(delta);
    
    % Store results
    std_d        = [std_d; stdVals];
    mean_d       = [mean_d; meanVals];
    mean_ab_d       = [mean_ab_d; [mean(test_results(:,1)), mean(test_results(:,2))]];
end

% Plot first set
figure;   % Make a new figure
errorbar(dims, mean_d, std_d, '.', 'MarkerSize', 20);
hold on;  % Keep the same figure for second plot

% % --- Second set (N=4) ---
% std_d2 = [];
% mean_d2 = [];
% sigma_d2 = [];
% mean_delta_d2 = [];
% dims2 = 4:6;
% 
% for d = dims2
%     filename = sprintf('random_state_d%d_N4_r%d.mat', d, d-1);
%     load(filename, 'test_results');
%     N = size(test_results, 1);
% 
%     % Compute mean
%     meanVals = sum(test_results, 1) / N;
%     mean_delta = meanVals(2) - meanVals(1);
% 
%     % Compute std and sigma
%     stdVals = [0 0];
%     sigma = 0;
%     for i = 1:N
%         stdVals(1) = stdVals(1) + (meanVals(1) - test_results(i,1))^2;
%         stdVals(2) = stdVals(2) + (meanVals(2) - test_results(i,2))^2;
%         sigma = sigma + (mean_delta - (test_results(i,2) - test_results(i,1)))^2;
%     end
% 
%     % Store results
%     std_d2        = [std_d2; stdVals];
%     mean_d2       = [mean_d2; meanVals];
%     sigma_d2      = [sigma_d2; sigma];
%     mean_delta_d2 = [mean_delta_d2; mean_delta];
% end
% 
% % Plot second set on top of the first
% errorbar(dims2, mean_delta_d2, sigma_d2, '.', 'MarkerSize', 20);

% Adjust axes, labels, etc.
xticks(3:10);
xlim([3.5 9.5]);
ax = gca;
ax.FontSize = 14;
xlabel('d','Interpreter','latex','FontSize',18)
ylabel('$\Delta$','Interpreter','latex','FontSize',18)

%legend('N=3','N=4','Location','Best');

hold off;


