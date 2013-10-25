function rgb = pbBitStealing2RGB(bsnumber,rgbMatrix,checkRange)
%
% function bsnumber = pbBitStealing2RGB(bsnumber,rgbMatrix,checkRange)

%
% Converts indices into RGB values using rgbMatrix.
% bsnumber   :   vector of bit-stealing indices
% rgbMatrix :   n x 3 matrix of rgb values
% checkRange :  if 1, will check bsindex to ensure that all values are in
%               range. Set to zero if you want to skip this step. Default is 1.
%

% 6-Oct-08	: altered code to return first 3 columns of rgbMatrix, rather than all columns - PJB

if nargin<3
    checkRange=1;
end;
if nargin<2
    error('insufficient inputs to function');
end;

[nr,nc]=size(bsnumber);
bsnumber=reshape(bsnumber,nr*nc,1);

if checkRange==1
    nrgb = size(rgbMatrix,1);
    bsLow = find(bsnumber<1);
    bsnumber(bsLow) = 1;
    bsHigh = find(bsnumber>nrgb);
    bsnumber(bsHigh)=nrgb;
end;


rgb=rgbMatrix(bsnumber,:);
rgb=reshape(rgb,nr,nc,3);
