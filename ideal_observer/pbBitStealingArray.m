function rgbValues = pbBitStealingArray(rgbIndex)
% function rgbValues = pbBitStealingArray
%
% makes a matrix of rgb combinations 
%
global pbBigRGBArray;

res = 7; % the number of stolen samples between points

numpts = 256;

% the code for incrementing r,g, and b
rinc = repmat([0 0 1 1 0 0 1],1,numpts-1)';
ginc = repmat([0 0 0 0 1 1 1 ],1,numpts-1)';
binc = repmat([0 1 0 1 0 1 0],1,numpts-1)';

% makes a list of integers from 1 to [numpts], repeating each for [res] times
basepix = ([1:numpts-1]'*ones(1,res))';
basepix = basepix(:);

% the values of r,g, and b to plug into pellipower
rpix = [basepix+rinc;numpts];
gpix = [basepix+ginc;numpts];
bpix = [basepix+binc;numpts];

% final data matrix, in columns
pbBigRGBArray = [rpix,gpix,bpix]-1;

if(exist('rgbIndex','var')==1)
	rgbValues=pbBigRGBArray(rgbIndex,:);
else
	rgbValues=pbBigRGBArray;
end
	
return
