%% This script produces figures of the stimuli for use in the paper.
% TSAW

clear all; close all;


%% Parameters used in the experiment:
wNZ = .005; 

% filtering params:
sfCF = [NaN 3.2]; % the cutoff frequencies for blur kernal.
sfBW = [NaN NaN NaN]; % not applicable.
% orientation CO
orientCF = [NaN NaN NaN];
% half-width BW
orientBW = [NaN NaN NaN];

o1 = orientCF - orientBW; o2 = orientCF + orientBW;

% load images used in the experiment:
load('~/Dropbox/Projects/AMD_Faces/clinical_faces_analysis/images/faceFiltStimuliOrderedViews.mat')

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

%% apply blur kernels to one face:
ident = 1;
viewPoints = [1:5,1:5];
viewpoint = viewPoints(ident);
tmp=find(personID==ident); % select subset of items that are the right person
curfaceindex = tmp(viewpoint);

testFace = faceViews(:,:,curfaceindex);

bigIm = NaN(nr*2,nc);

for blur = 1 : length(sfCF)
    fv=bpimage(testFace,theFilters(:,:,blur),0);
    
    startRow = (blur-1)*nr + 1;
    endRow = startRow + nr - 1;
    startCol = 1;
    endCol = nc;
    bigIm(startRow:endRow,startCol:endCol) = fv;
    
end
% figure;
% sc(bigIm);

saveIm = normVector(bigIm);
imwrite(saveIm,'~/Dropbox/Projects/AMD_Faces/clinical_faces_analysis/images/blur_cutoff.jpg','jpg','Quality',100)

