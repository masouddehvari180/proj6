% program for create .csv files which contaims both model and corresponding
% ERA5 data. the usage has been illustrated in 4 steps as below:

clear all
close all


% step 1:
addpath( genpath('E:\MATLAB\R2021b\toolbox\nctoolbox') ) % replace the folder of nctoolbox
setup_nctoolbox



% step 2:
str='J:\Grib2Files_202211'; % folder of *.grib files data

%step 3:
out_fol='H:\proj6\remained_cycles'; % determine output foder

a1=dir([str '/2022*']);

% step 4: 
nc_n='All2022_2023.nc'; % ERA5 data path and name



%----------------- main code -------------%

time=ncread(nc_n,'time');
T_era = datetime(1900,1,1) + hours(time);
T_era.Format = 'yyyy-MM-dd HH:mm:ss';
mjd2=mjuliandate(datevec(T_era));

lon_era=ncread(nc_n,'longitude');
lat_era=ncread(nc_n,'latitude');
u_era=ncread(nc_n,'u');
v_era=ncread(nc_n,'v');


lat_c  = 30:0.25:46;
lon_c  = 48:0.25:56 ;

[lon_c,lat_c]=meshgrid( lon_c,lat_c );
lon_c=lon_c(:);
lat_c=lat_c(:);


for ii=1:size(a1,1)
    
    
    data_fol=[a1(ii).folder '/' a1(ii).name];
    
    
    a=dir([ data_fol '/*.grb2' ]);
    if size(a,1)>=1
        
        
        BB1=zeros(1e6,11);
        i2=1;
        
        for i=1:size(a,1)
            
            size(a,1)-i
            f1=a(i).name;
            %
            t1=f1(5:15);
            t1=datetime( str2num( [t1(1:4) ' ' t1(5:6) ' ' t1(7:8) ' ' t1(9:10) ' 0 0' ]));
            t1.Format = 'yyyy-MM-dd HH:mm:ss';
            mjd1=mjuliandate(datevec(t1));
            [t_index,~]=find(mjd2==mjd1);
            name=f1(19:end-5);
            if isempty(t_index)
                continue;
            end
            f=[data_fol '/' f1];
            
            ds2 = ncdataset(f);
            longitude=ds2.data('lon');
            latitude=ds2.data('lat');
            % [longitude,latitude]= meshgrid(longitude,latitude);
            isobaric=ds2.data('pressure');
            u2 = double(ds2.data('U-component_of_wind'));
            u2=squeeze(u2);
            v2 = double(ds2.data('V-component_of_wind'));
            v2=squeeze(v2);
            T= double( ds2.data('Temperature') );
            T=squeeze(T);
            q= double( ds2.data('Specific_humidity') );
            q=squeeze(q);
            msl_p= double( ds2.data('Pressure_msl') );
            msl_p=squeeze(msl_p);
            
            [o1,~]=find( isobaric==100000 );
            if isempty(o1)
                error('isobaric level could not found.')
            end
            
            u2=squeeze( u2(o1,:,:) );
            v2=squeeze( v2(o1,:,:) );
            T=squeeze( T(o1,:,:) );
            q=squeeze( q(o1,:,:) );
            
            %-------------- interpolation in grib2 -------------%
            
            
            u_f= interp2( longitude,latitude,u2,lon_c,lat_c );
            v_f= interp2( longitude,latitude,v2,lon_c,lat_c );
            T_f= interp2( longitude,latitude,T,lon_c,lat_c );
            q_f= interp2( longitude,latitude,q,lon_c,lat_c );
            msl_f= interp2( longitude,latitude,msl_p,lon_c,lat_c );
            
            %------------- interp era ---------------%
            
            u_f_era= interp2( lat_era , lon_era , u_era(:,:,t_index) , lat_c , lon_c );
            v_f_era= interp2( lat_era , lon_era , v_era(:,:,t_index) , lat_c , lon_c );
            
            %-------------------------------%
            
            m1=length(lat_c);
            index=[1:m1]';
            time1=repmat(t1,[m1,1]);
            %BB  = cell2table(cell(0,11), 'VariableNames', {'index', 'Lon', 'Lat', 'Time', 'Msl_P', 'T', 'q', 'U', 'V','U_era','V_era'});
            
            %     BB{1,1}(i2:i2+m1-1,1)=index;
            %     BB{1,2}(i2:i2+m1-1,1)=lon_c;
            %     BB{1,3}(i2:i2+m1-1,1)=lat_c;
            BB{1,4}(i2:i2+m1-1,1)=(time1);
            %     BB{1,5}(i2:i2+m1-1,1)=msl_f;
            %     BB{1,6}(i2:i2+m1-1,1)=T_f;
            %     BB{1,7}(i2:i2+m1-1,1)=q_f;
            %     BB{1,8}(i2:i2+m1-1,1)=u_f;
            %     BB{1,9}(i2:i2+m1-1,1)=v_f;
            %     BB{1,10}(i2:i2+m1-1,1)=u_f_era;
            %     BB{1,11}(i2:i2+m1-1,1)=v_f_era;
            BB1(i2:i2+m1-1,1)=index;
            BB1(i2:i2+m1-1,2)=lon_c;
            BB1(i2:i2+m1-1,3)=lat_c;
            %     BB1(i2:i2+m1-1,4)=(time1);
            BB1(i2:i2+m1-1,5)=msl_f;
            BB1(i2:i2+m1-1,6)=T_f;
            BB1(i2:i2+m1-1,7)=q_f;
            BB1(i2:i2+m1-1,8)=u_f;
            BB1(i2:i2+m1-1,9)=v_f;
            BB1(i2:i2+m1-1,10)=u_f_era;
            BB1(i2:i2+m1-1,11)=v_f_era;
            
            
            
            
            
            
            i2=i2+m1;
        end
        BB1(i2:end,:)=[];
        
        
        
        
        
        
        
        
        T={'index', 'Lon', 'Lat', 'Time', 'Msl_P', 'T', 'q', 'U', 'V','U_era','V_era'};
        m1=size(BB1,1);
        
        for i=1:size(T,2)
            if i==4
                T(2:m1+1,i)=num2cell(BB{1,i});
            else
                T(2:m1+1,i)=num2cell(BB1(:,i));
            end
            
        end
        
        out_name1=[out_fol '/Cycle_' name '.csv'];
        writecell(T,out_name1)
        T2= sort_tbl1(out_name1);
        writetable(T2,out_name1)
        
        clear BB T T2
    end
    
end



function T2= sort_tbl1(f)



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




end









































