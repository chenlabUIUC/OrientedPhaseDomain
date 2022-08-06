function [fitresult, gof] = createFit(pos_x, pos_y)
%CREATEFIT(POS_X,POS_Y)
%  Create a fit.
%
%  Data for 'untitled fit 1' fit:
%      X Input : pos_x
%      Y Output: pos_y
%  Output:
%      fitresult : a fit object representing the fit.
%      gof : structure with goodness-of fit info.
%
%  See also FIT, CFIT, SFIT.

%  Auto-generated by MATLAB on 06-Aug-2022 02:39:43


%% Fit: 'untitled fit 1'.
[xData, yData] = prepareCurveData( pos_x, pos_y );

% Set up fittype and options.
ft = fittype( 'gauss4' );
opts = fitoptions( 'Method', 'NonlinearLeastSquares' );
opts.Display = 'Off';
opts.Lower = [0 0.61 0 0 0.65 0 0 0.68 0 0 0.71 0];
opts.StartPoint = [427 0.6925 0.000746119263031986 383.999959332778 0.6955 0.000723724752647086 167.999959332778 0.6895 0.00150381155494339 110.999986759439 0.6985 0.00105535315618804];
opts.Upper = [0 0.65 Inf 0 0.68 Inf Inf 0.71 Inf 0 0.76 Inf];

% Fit model to data.
[fitresult, gof] = fit( xData, yData, ft, opts );

% Plot fit with data.
figure( 'Name', 'untitled fit 1' );
h = plot( fitresult, xData, yData );
legend( h, 'pos_y vs. pos_x', 'untitled fit 1', 'Location', 'NorthEast' );
% Label axes
xlabel pos_x
ylabel pos_y
grid on


