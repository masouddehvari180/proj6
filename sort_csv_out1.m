clear all
close all
%
data_fol='J:\cycle_out';
out_fol= 'H:\proj6\converted_csv_files';

a=dir([data_fol '\*.csv']);

for ii=1:size(a,1)
    
    size(a,1)-ii
    f1= a(ii).name;
    
    
    
    
    f=[data_fol '\' f1];
    
    T = readtable(f);
    %
    
    mjd1=mjuliandate( datevec(T.Time) );
    mjd1= mjd1-mjd1(1);
    T.TfFd= mjd1;
    T= [T(:,1:4)  T(:,12) T(:,5:11)];
    
    %
    time=T.Time;
    index= T.index;
    lat= T.Lat;
    lon= T.Lon;
    msl= T.Msl_P;
    t= T.T;
    q= T.q;
    u= T.U;
    v= T.V;
    %------------ sort by index --------%
    
    T2=T(1,:);
    [m1,m2]=size(T);
    % T2(m1,m2)=0;
    T2=[];
    ind1=unique(index);
    
    for i=1:length(ind1)
        
        [o1,~]=find(index==ind1(i));
        
        qw=T(o1,:);
        
        mjd1=mjuliandate(datevec(qw.Time));
        [o1,o2]=sort(mjd1,'ascend');
        qw=qw(o2,:);
        
        T2=[T2;qw];
        
    end
    
    writetable(T2,[out_fol '/' f1])
   
    
end



