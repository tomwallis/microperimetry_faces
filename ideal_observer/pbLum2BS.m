function bitsValues=pbLum2BS(lumValues,L,B)

    if(nargin<3)
        error('must pass 3 variables to this function');
    end;
    if(size(L,2)~=1)
        L=reshape(L,prod(size(L)),1);
    end;
    if(size(B,2)~=1)
        B=reshape(B,prod(size(B)),1);
    end;
    if(length(L)~=length(B))
        error('L and B must have the same number of elements');
    end;
    [nr,nc]=size(lumValues);
    lv1=reshape(lumValues,nr*nc,1);
    Lmin = min(L);
    Lmax = max(L);
    tmpMin = find(lv1<Lmin);
    tmpMax = find(lv1>Lmax);
    if (length(tmpMin)>0)
        lv1(tmpMin) = Lmin;
    end;
    if (length(tmpMax)>0)
        lv1(tmpMax) = Lmax;
    end;
   
    bitsValues=interp1(L,B,lv1,'nearest');
%     bitsValues=interp1q(L,B,lv1);
    bitsValues=reshape(bitsValues,nr,nc);
    
    
