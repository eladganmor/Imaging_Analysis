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
[fluo, r, rate, stim] = genCaTrace(stimMat,tuning,bias);%generate fluorescence trace, spike counts, and rate

%% Fit

data = calcFValsSimple(fluo,stim);
params = searchBestParamsSimple(data);
tuningFit = exp(params(1) + params(2:end));

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