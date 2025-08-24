power = load("power.mat").Power;

% check if the data is normally distributed using histograms
subplot(1,2,1)
histogram(power(1, :))
title('Successful Condition')

subplot(1,2,2)
histogram(power(2, :))
title('Unsuccessful Condition')

% t-test to compare mean
[h, p] = ttest2(power(1, :), power(2, :));

if h == 0
    disp('The means are not significantly different.');
else
    disp('The means are significantly different.');
end