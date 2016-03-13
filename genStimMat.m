nCycles = 3e3;%how many overall stimulus presentations
nStim = 9;

tOn = 0.333;%on duration of stimulus
tOff = 0.167;%inter stimulus interval duration of stimulus


stimLength = round(tOn/presets.dt);
isi = round(tOff/presets.dt);
T = nCycles*(isi+stimLength) + isi;
stimStart = isi+1 : isi + stimLength : T - isi;
stimMat = zeros(nStim,T);
for cycleInd = 1:nCycles
    s = randi(nStim);
    if rand>0.1
        stimMat(s, stimStart(cycleInd) : stimStart(cycleInd) + stimLength - 1) = 1;
    end
end
stimMat(:,1:100) = 0;%start the experiment with a blank

%%Stimulus tuning
k = 0.1;
mu = pi/2;
b = -2;
stimVals = linspace(0,pi,nStim);
tuning = exp(logVonMises(stimVals,[b,mu,k]));
