load('Y:\nick\behavior\grooming\1p\outputs\long_groom_global_5spre.mat')
fs = 90;
blen = 5;
clear yFitted

all_global = all_global(~cellfun('isempty',all_global));
clc
modelfun = @(b,x) b(1) * exp(-b(2)*x(:, 1));

% close all

for i = 1:length(all_global)
    X = all_global{i}(blen*90+1:end);
    t = xt(X, fs);

    tbl = table(t', X');

    beta = [max(X) 0.5];

    try
        mdl = fitnlm(tbl, modelfun, beta);
        R_squared(i) = mdl.Rsquared.Ordinary;

        coeffs(:,i) = mdl.Coefficients.Estimate;

        coefficients = mdl.Coefficients{:, 'Estimate'};
        yFitted{i} = coefficients(1) * exp(-coefficients(2)*t);
        % Compute the residuals
        residuals = X - yFitted{i};

        
        
%         % Compute the total sum of squares (SS_tot)
%         SS_tot = sum((X - mean(X)).^2);
%         
%         % Compute the sum of squared residuals (SS_res)
%         SS_res = sum(residuals.^2);
%         
%         % Compute R-squared
%         R_squared(i) = 1 - (SS_res / SS_tot);

        figure, plot(t, X);
        hold on, plot(t, yFitted{i})
        title(['a = ', num2str(coeffs(1,i)),';   b = -',num2str(coeffs(2,i)), 'R-squared: ', num2str(R_squared(i))])
    catch
        disp('Could not fit an exponential decay - Fit a line instead')
        p{i} = polyfit(t, X, 1);
    end


    
end