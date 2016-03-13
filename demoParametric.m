%% set presets
global presets

presets.dt = 1/6;
presets.stimHistoryLength = presets.dt;
%%parameters governing trasformation from spikes to calcium
presets.tau = 0.25;
presets.a = 1;
presets.rMax = 10;
presets.sig = 0.5;
presets.rMax = 20;%if we underestimate sigma, we need this

genStimMat;
bias = -1;
%% run simulation
decayFun = presets.a*exp(-(0:100)*presets.dt/presets.tau);
tuningFun = @logVonMises;
initParams = [0,0,1];%baseline, mu, k

[fluo, r, rate, stim] = genCaTrace(stimMat,tuning,bias);%generate fluorescence trace, spike counts, and rate
theta = stimVals*stim;%this just means that we get the value of the stimulus on the screen at each time
theta(sum(stim)==0) = NaN;%blank

%% Fit

data = calcFValsSimple(fluo,theta);
data.xVals = data.fVals(2,:);%remove bias feature, it should be taken care of by tuning function
[params,~,rateEst] = searchBestParamsSimple(data,initParams,tuningFun);
tuningFit = exp(tuningFun(stimVals,params));

%% plot
figure
subplot(1,3,1:2)
t = (1:length(fluo))*presets.dt;
plot(t,fluo,'k')
xlabel('Time (s)')
ylabel('Fluorescence')

subplot(1,3,3)
plot(stimVals*180/pi, exp(bias+tuning)/presets.dt,'k')
hold on
plot(stimVals*180/pi, tuningFit/presets.dt,'b')
xlabel('Orientation (deg)')
ylabel('Rate (spikes/s)')
legend('True','Fit')