%we have to figure out how the range you choose to fit would effect your age result; another issue is how to determine the temperature

clear all;
close all;
%========time conversion===================================================
time2Year = 1/(10^12*3600*24*365);

%% Plot the Barium and Strontium concentration vs. distance, with data point

% Load data from .mat file, then we can get arrays Sr, Ba, microns
temp = 820 + 273;
s = load('ESS134_Data_1.mat','-mat');
microns = s.microns;

Ba = s.Ba(:,1);
BaError = s.Ba(:,2);

Sr = s.Sr(:,1);
SrError = s.Sr(:,2);

figure;
plot(microns, Sr, '--rs', 'LineWidth', 1, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize', 10 );
title('Sr  Concentration vs. x');
xlabel('x(microns)');
ylabel('concentration(ppm)');

figure;
plot(microns, Ba, '--rs', 'LineWidth', 1, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize', 10 );
title('Ba  Concentration Fitting Curve');
xlabel('x(microns)');
ylabel('concentration(ppm)');
hold on;



%% Using non-linear fit to fit the Ba experimental data

% Function model: a(1) = 0.5*(c_left+c_right); a(2) = 0.5*(c_right-c_left); a(3) = D*t;
% concentration = @(a, X) 1000*a(1) + 1000*a(2)*erf((X(:,1)-a(3))/(2*sqrt(a(4))));


%=====================Using nlinfit to fit the data========================
concentration = @(a, X) 1000*a(1) + 1000*a(2)*erf((X(:,1)-a(3))/(2*a(4)));
BaBeta0 = [5, 1.5, 25, 5];

obs = [microns];

opts = statset('MaxIter', 10000);

opts.RobustWgtFun = 'bisquare';
% fitBegin = 60;------to reproduce 1073 years
% fitEnd = 120;------ to reproduce 1073 years
fitBegin = 76;
fitEnd = 90;
%how the fitting affects time
%result = lsqcurvefit(concentration, beta0,  microns(60:100, 1), Ba(60:100, 1));
%res        ult = nlinfit(obs(60:100, 1), Ba(60:100, 1), concentration, beta0, opts);


%===================Calculation Of uncertainty=============================
[BaResult, resid, J, sigma] = nlinfit(obs(fitBegin:fitEnd, 1), Ba(fitBegin:fitEnd, 1), concentration, BaBeta0, opts);
ci = nlparci(BaResult, resid, 'covar', sigma);


%===================Plot the fit data======================================
%BaFit = 1000*result(1) + 1000*result(2)*erf((microns(60:100, 1)-result(3))/(2*sqrt(result(4))));
%BaFit = 1000*BaResult(1) + 1000*BaResult(2)*erf((microns(fitBegin:fitEnd, 1)-BaResult(3))/(2*sqrt(BaResult(4))));
BaFit = 1000*BaResult(1) + 1000*BaResult(2)*erf((microns(fitBegin:fitEnd, 1)-BaResult(3))/(2*BaResult(4)));
%plot(microns(60:100, 1), BaFit,  'k', 'LineWidth', 2);
plot(microns(fitBegin:fitEnd, 1), BaFit,  'k', 'LineWidth', 2);


%===============Plot the inversion of erf, straight line fit===============
figure;

invBa = erfinv((Ba - 1000 * BaResult(1)) / (1000 * BaResult(2)));
invBaFit = erfinv((BaFit - 1000 * BaResult(1)) / (1000 * BaResult(2)));

plot(microns(fitBegin:fitEnd), invBa(fitBegin:fitEnd), '--rs', 'LineWidth', 1, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize', 10 );
title('Inverse ERF Linear Fit');
xlabel('x(microns)');
ylabel('erfinv((Concentration-C_1)/C_2)');
hold on;
plot(microns(fitBegin:fitEnd), invBaFit, 'k', 'LineWidth', 2);


%===================Why I select this range for nonlinear fitting================================
figure;
invBa = erfinv((Ba - 1000 * BaResult(1)) / (1000 * BaResult(2)));
plot(microns(:), invBa, '--rs', 'LineWidth', 1, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize', 10);
title('Inverse ERF Linear Fit(Large Data Range)');
xlabel('x(microns)');
ylabel('erfinv((Concentration-C_1)/C_2)');


%==================Goodness of fit using chi-square value---------------------
chisquareBa = sum( ( (Ba(fitBegin:fitEnd)' - BaFit(:)').^2)./ (BaError(fitBegin:fitEnd)'.^2) ) / (15-1-4); 


%===================Age Calculation----------------------------------------------
diffusivityBa = 2.9 * 0.1 * exp(-455000 / (8.314472 * 1093));
% diffusivityBa = 2.9 * 0.1 * exp(-450000 / (8.3 * 1073));-------reproduce
% 1073 yers

ageBa = BaResult(4)^2 / diffusivityBa * 10^-12 / (3600 * 24 * 365);
disp(ageBa);


%==================Age Uncertainty Calculation(2-sigma)-----------------------------
% DTSqrtBa = BaResult(4);
% DTSqrtBaUnty = (ci(4,2)-ci(4,1))/2;
% gasConstant = 8.314472;
% diffusivityBaUnty = sqrt((diffusivityBa*20000 / (gasConstant*temp))^2 + (diffusivityBa*455000*20 / (gasConstant*temp^2))^2);
% 
% ageBaUnty = sqrt((DTSqrtBa/diffusivityBa^2 * diffusivityBaUnty)^2 + (DTSqrtBa/diffusivityBa * DTBaUnty)^2);
% disp(ageBaUnty/365/24/3600);
%=================age uncertainty recalculation----------------------------
% Dmin = 2.9 * 0.1 * exp(-475000 / (8.314472 * 1073));
% Dmax = 2.9 * 0.1 * exp(-435000 / (8.314472 * 1113));
% % DTmin = ci(4,1);
% % DTmax = ci(4,2);
% ageMax = DTSqrtBa / Dmin / 10^12 / (3600*24*365);
% ageMin = DTSqrtBa / Dmax / 10^12 / (3600*24*365);
%=================Monte Carlo method to calculate age uncertainty---------------------------------------
nRandom = 10^5;
BaResult4Unty = ci(4,2)-BaResult(4);
nBin = 200;

numeratorRan = normrnd(455000, 20, nRandom, 1);
tempRan = normrnd(1093, 10, nRandom, 1);
DTSqrtRan = normrnd(BaResult(4), BaResult4Unty, nRandom, 1);
% DTSquareRan = DTSqrtRan.*DTSqrtRan;
diffusivityRan = 2.9 * 0.1 * exp(-numeratorRan./ (8.314472 * tempRan));
% ageRan = sqrt(DTSquareRan./(diffusivityRan.^2));
ageRan = DTSqrtRan.^2./diffusivityRan*time2Year;

figure;
hist(ageRan, nBin);
title('Age Monte Carlo Histogram');
xlabel('age(year)');
ylabel('number of events');

histHandle = hist(ageRan, nBin);
yearsPerBin = max(ageRan)/nBin;
ageRanMean = mean(ageRan);
%===============integral histogram-----------------------------------------
cumAgeRan = cumsum(histHandle)/sum(histHandle);
figure;
plot(yearsPerBin*[0:nBin], [0, cumAgeRan]);
title('Age CDF');
xlabel('age(year)');
ylabel('probablity');

ageSum = sum (ageRan);
hold on;
scatter(155.76, 0.05, 'fill');
hold on;
scatter(2959.6, 0.95, 'fill');

%% -----------------Calculation of Sr diffusivity------------------------
concentrationSr = @(a, X) a(1) + a(2)*erf((X(:,1)-a(3))/(2*sqrt(a(4))));

SrBeta0 = [80, 20, 25, 0.8];

obs = [microns];

opts = statset('MaxIter', 10000); 
opts.RobustWgtFun = 'bisquare';

fitBeginSr = 70;
fitEndSr = 100;

SrResult = nlinfit(obs(fitBeginSr:fitEndSr, 1), Sr(fitBeginSr:fitEndSr, 1), concentrationSr, SrBeta0, opts);
SrFit = SrResult(1) + SrResult(2)*erf((microns(fitBeginSr:fitEndSr, 1)-SrResult(3))/(2*sqrt(SrResult(4))));

figure;
plot(microns, Sr, '--rs', 'LineWidth', 1, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize', 10 );
hold on;
plot(microns(fitBeginSr:fitEndSr, 1), SrFit,  'k', 'LineWidth', 2);

SrFit = SrResult(1) + SrResult(2)*erf((microns(fitBegin:fitEnd, 1)-SrResult(3))/(2*sqrt(SrResult(4))));

figure;
plot(microns, Sr, '--rs', 'LineWidth', 1, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'g', 'MarkerSize', 10 );
title('Sr  Concentration Fitting Curve');
xlabel('x(microns)');
ylabel('concentration(ppm)');
hold on;
plot(microns(fitBegin:fitEnd, 1), SrFit,  'k', 'LineWidth', 2);

%=============Recalculate the diffusivity of Sr---------------------------
SrDT = SrResult(4);
diffusivitySr = SrDT / ageRanMean;

%% Calculate the goodness of fit using the Chi-Squre misfit Value
%freedom = numel(microns) - 4 -1;
%reducedBaChiSquare = (BaFit - Ba)^2./(BaError * freedom);

%reducedBaChiSquare = (SrFit - Sr)^2./(SrError * freedom);