function rgb = pbBSindex2RGB(bsindex)
%
% function bsnumber = pbBSindex2RGB(bsindex)
%
% Converts an index into an RGB value using the Bit-stealing array
%
global pbBigRGBArray;

if (size(bsindex,2)~=1)
	error('bsindex input must be a vector');
end;

if isempty(pbBigRGBArray)
	tmp=pbBitStealingArray;
else
	tmp=pbBigRGBArray;
end;

tmpMax=max(bsindex(:));
tmpMin=min(bsindex(:));

if tmpMin<1
    error('minimum index cannot be less than 1');
end;
if tmpMax>pbMaxBITStealingIndex
    error('maximum bit-stealing index is out of bounds');
end;

n=length(bsindex);
rgb=zeros(n,3); % assumes RGB is a 3 colum matrix
for kk=1:n
	rgb(kk,:)=tmp(bsindex(kk),1:3);
end;



