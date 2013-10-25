%% This script implements an ideal observer for the 10AFC faces experiment.
% Will return performance of the ideal observer for a number of simulated
% trials at each blur level (percent correct). Duplicated for 3 "blocks".
%
% The ideal observer responds with the face that has the largest
% cross-correlation with the target (see Gold et al 1999).
%
% TSAW

clear all; close all;

appendData = 0;  % if 1, append new trials to existing data file.


%% Parameters for the ideal observer analysis:
nTrialsPerBlurLevel = 5; % the number of trials per condition to have the observer do.

% here, condition = 3 blur levels x N levels of white noise contrast. Randomly select faces and
% views.

% set random stream:
rand("state","reset")


%% Parameters used in the experiment:
wNZ = [0.005 logspace(log10(.5),log10(1000),13)]; % the variance on the Gaussian pixel noise.


if appendData == 1 && exist('idealPerformance.mat','file')==2
    load('idealPerformance.mat')
else
    dataMat = [];
end


% filtering params:
sfCF = [NaN 32 16]; % the cutoff frequencies for blur kernal.
sfBW = [NaN NaN NaN]; % not applicable.
% orientation CO
orientCF = [NaN NaN NaN];
% half-width BW
orientBW = [NaN NaN NaN];

o1 = orientCF - orientBW; o2 = orientCF + orientBW;

% load images used in the experiment:
load('/Users/tomwallis/Dropbox/Projects/AMD_Faces/clinical_faces_analysis/images/faceFiltStimuliOrderedViews.mat')

% set up image filters for blur kernals as used in the experiment (mostly taken from Chris' code faceSum.m):
% orientation filter parameters:
numFilters = length(sfCF);
[nr, nc] = size(faceViews(:,:,1)); % get size of the images
theFilters = zeros(nr,nc,numFilters);
for kk=1:numFilters
    [sf1(kk), sf2(kk)] = calcCutOff(sfCF(kk),sfBW(kk));
    theFilters(:,:,kk)=filter2d(0,sfCF(kk),o1(kk),o2(kk),nr);
end;

noiseContrast = wNZ;

% to make the images as similar to the experiment as possible, load all
% calibration data so we convert to luminance:
scrinfo.calfile ='eeg_1024x768x100Hz_Nov10.mat';
disp('reading calibration file...');
calfitrec=pbReadCalibrationFile('eeg_1024x768x100Hz_Nov10.mat');
rgbMat = calfitrec.caldata(:,1:3); % rgb matrix
B = calfitrec.caldata(:,4); % array of bit-stealing numbers
L = calfitrec.caldata(:,5); % array of luminance values

% EDIT: If the demo is to be run on a machine w/two displays change to 1.
%
% number = 0; for testing, 1 for experiment
scrinfo.number=0;
scrinfo.width=calfitrec.displaypixels(1);
scrinfo.height=calfitrec.displaypixels(2);
scrinfo.framerate=calfitrec.framerate;
scrinfo.pixelsize=calfitrec.pixelsize;
scrinfo.pixpercm=mean(calfitrec.displaypixels ./ calfitrec.displaycm);
dcNumber = calfitrec.backgroundNumber;
RGBgrey=pbBSindex2RGB(dcNumber);
avgLum = pbBitStealing2Lum(calfitrec.backgroundNumber,L,B);

%% Set up a trial structure, for three blocks at each blur level. I
% understand this could be done more efficiently without loops...

%matlabpool open

tic;

for noiseContrastLevel = 1 : length(noiseContrast)
    
    for blur = 1 : 3
        
        % increment data mat per blur level to allow parallelisation.
        % new rows:
        dataMatNew = NaN(nTrialsPerBlurLevel,6); % cols = blur level, wNZ, normalisedNoise, face ident, correct / incorrect.
        startRow = size(dataMat,1) + 1;
        endRow = startRow + nTrialsPerBlurLevel - 1;
        dataMat = [dataMat; dataMatNew];
        clear('dataMatNew');
        
        trialRowIndex = startRow:endRow;
        
        tmpDat = [];
        
        for trial = 1 : nTrialsPerBlurLevel
            
            % Select a face identity (1--10):
            curPerson = randi(10);
            
            % Select a view by selecting a subset of items that are the
            % current person:
            tmp=find(personID==curPerson); % select subset of items that are the right person
            rf = randi(length(tmp)); % select a random member of the the subset
            
            curFaceIndex = tmp(rf);
            curCondition = blur;
            curContrast = .05;
            curNoiseContrast = noiseContrast(noiseContrastLevel);
            
            %% select image; apply filtering and noise:
            fv=bpimage(faceViews(:,:,curFaceIndex),theFilters(:,:,curCondition),0);
            curFaceImage =  fv.*(curContrast)./sqrt(var(fv(:)));  % extract image & set to correct contrast
            
            % Set up nr x nc noise
            if(curNoiseContrast > 0)
                %curNoise = bpimage(randn(nr),theFilters(:,:,curCondition),0);
                curNoise = randn(nr);
                sdtmp = sqrt(var(curNoise(:)));
                curNoise = curNoise .* (curNoiseContrast/sdtmp);
                upperLim=2*sqrt(curNoiseContrast);
                lowerLim=-upperLim;
                % remove out-of-range values from curNoise
                tmp=find((curNoise > upperLim)|(curNoise<lowerLim));
                if(~isempty(tmp))
                    newNz=sqrt(curNoiseContrast)*randn(length(tmp),1);
                    outOfRange=find((newNz > upperLim)|(newNz<lowerLim));
                    while(~isempty(outOfRange))
                        newNz(outOfRange)=sqrt(curNoiseContrast)*randn(length(outOfRange),1);
                        outOfRange=find((newNz > upperLim)|(newNz<lowerLim));
                    end;
                    curNoise(tmp)=newNz;
                end;
            else
                curNoise = zeros(nr,nc);
            end;
            
            %Create Face+Noise contrast array
            stimPlusNoise=curFaceImage+curNoise;
            
            % convert from contrast to luminance:
            stimLum=avgLum*(1+stimPlusNoise);
            
            % convert Luminance to BitStealing Numbers:
            stimBSnumbers=pbLum2BS(stimLum,L,B);
            
            % finally, convert BS numbers to RGB values:
            stimRGBarray = pbBitStealing2RGB(stimBSnumbers,rgbMat,0);
            
            %% calculate the ideal observer's response:
            
            % retrieve the 10 response alternatives for this viewpoint:
            i = 1:10;
            responseImNums = (i-1)*5 + rf;
            
            crossCorr = zeros(1,10);
            
            for resp = 1 : 10
                %% select image; apply filtering and noise:
                rIm=bpimage(faceViews(:,:,responseImNums(resp)),theFilters(:,:,curCondition),0);
                rIm=rIm.*(curContrast)./sqrt(var(rIm(:)));  % extract image & set to correct contrast
                
                % convert to lum:
                rImLum=avgLum*(1+rIm);
                
                % convert Luminance to BitStealing Numbers:
                rImBSnumbers=pbLum2BS(rImLum,L,B);
                
                % finally, convert BS numbers to RGB values:
                rImRGBarray = pbBitStealing2RGB(rImBSnumbers,rgbMat,0);
                
                % calculate the cross correlation between the stimulus face and
                % the possible response face:
                %                 crossCorr(resp) = corr2(stimPlusNoise,respIm(:,:,resp));
                crossCorr(resp) = corr2(stimBSnumbers,rImBSnumbers);
            end;
            
            % the ideal observer responds with the maximum cross
            % correlation:
            [val, index] = max(crossCorr);
            
            % correct for possible tie by choosing one randomly:
            if length(index) > 1
                index = index(randi(length(index)));
            end
            
            % which face did observer pick?
            decision = responseImNums(index);
            
            % correct or incorrect?
            correct = decision==curFaceIndex
            
           
            tmpDat1 = [blur, curNoiseContrast, std(curNoise(:)), std(curFaceImage(:)), curPerson, correct];
            
            tmpDat(trial,:) = tmpDat1;

        end
        
        dataMat(trialRowIndex,:) = tmpDat;
    end
    
end

timeInHours = toc/60/60;

%matlabpool close

% sc(stimPlusNoise,[-1,1]);


save('idealPerformance_Octave.mat','dataMat')