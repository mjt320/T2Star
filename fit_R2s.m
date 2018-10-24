function [S0,R2s,T2s,RSq,model] = fit_R2s(TE_s,magnitude,isFit,threshold1,threshold2);
% Fit a single T2/T2* decay
% OUTPUT:
% S0: estimate of signal amplitude for TE=0
% R2s: estimate of R2/R2*
% T2s: estimate of T2/T2*
% RSq: R-squared value for fit
% model: vector of signals predicted by model
% INPUT:
% TE_s: vector of TE times
% magnitude: vector of magnitude signals
% isFIT: vector of ones and zeros determining which signals should be fit
% threshold1: only do fitting if first signal is greater than this (to prevent fitting of air voxels)
% threshold2: only include signals greater than this value in fitting  (to prevent fitting of noise at long TE)

tic;

S0=nan(1,1);
R2s=nan(1,1);
T2s=nan(1,1);
RSq=nan(1,1);
model=nan(size(magnitude));

s=@(coef,t) coef(1)*exp(-coef(2)*t); % this is the signal model function

if magnitude(1)<threshold1; return; end %skip low intensity voxels

t=TE_s(isFit==1); %define independent variable (time)
y=squeeze(magnitude(isFit==1)); %define dependent variable (signal)

%% use linear regression to determine initial guess
if size(y,1)>2; NLinReg=3; else NLinReg=2; end
temp=regress(log(y(1:NLinReg)),[ones(NLinReg,1) t(1:NLinReg)]);
x0=[exp(temp(1)) -temp(2)];
if x0(1)<0 || imag(x0(1))~=0 || isnan(x0(1)) || isinf(x0(1))...
        || x0(2)<0 || imag(x0(2))~=0 || isnan(x0(2)) || isinf(x0(2)); x0=[1000 10]; end


%% use non-linear least squares to fit
typicalX=[max(y) 50];

% use only data points above opts.threshold2 to avoid fitting noise
% if there are fewer than 2 data points selected then don't fit
idx_reduced=y>threshold2;
if sum(idx_reduced)<2; return; end

[x,resnorm,residual,exitflag,output]=lsqcurvefit(s,x0,t(idx_reduced),y(idx_reduced)...
    ,[],[],optimset('Display','off','TypicalX',typicalX,'Algorithm','levenberg-marquardt'));

S0=x(1); R2s=x(2);
model(idx_reduced)=s(x,t(idx_reduced));
RSq=1 - sum((y(idx_reduced).'-squeeze(model(idx_reduced)).').^2) / (sum((y(idx_reduced).'-mean(y(idx_reduced))).^2));

timeElapsed=toc;

if rand<0.001 %randomly plot data to check it's working
    figure(1),plot(t,y,'ko',t(idx_reduced),squeeze(model(idx_reduced)),'b-')
    title({['Initial coefficients: ' num2str(x0)] ['Fitted coefficients: ' num2str(x)] ['Time elapsed: ' num2str(timeElapsed)] ['Exit flag: ' num2str(exitflag)] ['Evaluations: ' num2str(output.funcCount)]});
    pause(0.1);
end

T2s=1/R2s;

end