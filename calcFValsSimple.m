%%input: 
%1) r - the fluorescence mearurement; 1xT
%2) stimMat - the stimulus matrix; nStimulixT

%%output:
%data - a struct containing the features for the optimization procedure

%this function relies on the global variable "presets"
function data = calcFValsSimple(r,stimMat)
    global presets
    data.response = single(r(1,:));%%In order to align data.fVals(t) with data.response(t+1)
    
    histLength = round(presets.stimHistoryLength/presets.dt);%length of stimulus history effects
    nStim = size(stimMat,1);
    data.fVals = zeros(nStim*histLength, size(stimMat,2), 'single');
    for hInd = 1:histLength%%create shifted copies of the stimulus (i.e. with delays)
        data.fVals((hInd-1)*nStim + 1 : hInd*nStim, :) = circshift(stimMat, hInd, 2);%%notice that even the first index is shifted by one relative to the stimulus, so we have no instantaneous effects
    end
    %%chop off the beginnig because of the circular shifting (the values
    %%there are not real)
    data.fVals = data.fVals(:,histLength+1 : end);
    data.response = data.response(histLength+1 : end);
    data.fVals = [ones(1,size(data.fVals,2)); data.fVals];%add a bias term
end