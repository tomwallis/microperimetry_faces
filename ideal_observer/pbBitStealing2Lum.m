function lumValues=pbBitStealing2Lum(bsValues,L,B,pellicoef)

    if (exist('pellicoef','var')==1)
        lumValues = pbPelliPower(pellicoef,bsValues);
    elseif (exist('L','var')==1)&(exist('B','var')==1)
        if(size(L,2)~=1)
            L=reshape(L,prod(size(L)),1);
        end;
        if(size(B,2)~=1)
            B=reshape(B,prod(size(B)),1);
        end;
%         L=reshape(L,1,prod(size(L)));
%         B=reshape(B,1,prod(size(B)));
        if(length(L)~=length(B))
            error('L and B must have the same number of elements');
        end;
        if(length(L)<2)
            lumValues=NaN;
            return;
        end;
        [nr,nc]=size(bsValues);
        bv1=reshape(bsValues,nr*nc,1);
        minB=min(B);
        maxB=max(B);
        tmp=find(bv1<minB); bv1(tmp)=minB;
        tmp=find(bv1>maxB); bv1(tmp)=maxB;
%         lumValues=interp1(B,L,bv1,'linear');
        lumValues=interp1q(B,L,bv1);
        lumValues=reshape(lumValues,nr,nc);
    else
        lumValues=NaN;
    end;
    
