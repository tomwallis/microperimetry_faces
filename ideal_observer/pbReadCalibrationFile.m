function calfitrec=pbReadCalibrationFile(path2file)
%
% function calfitrec=pbReadCalibrationFile(path2file)
%

    tmp=load(path2file);
    fn=fieldnames(tmp);
    calfitrec=getfield(tmp,fn{1});
end