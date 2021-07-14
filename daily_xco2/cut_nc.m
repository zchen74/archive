cd E:\co2_matlab
%Z.Chen, Oct 7th, cut the data into pieces based on time%
%Note this modified script only reads data north of 50 S and only take Land
%Nadir data, refer tp 
%filename='E:\co2_matlab\OCO2_b10c_10sec_GOOD_r3.nc4';

filename='OCO2_b10c_10sec_GOOD_r3.nc4';
si=ncread(filename,'sounding_id');
datt=ncread(filename,'date');
lonn=ncread(filename,'longitude');
latt=ncread(filename,'latitude');
qf=ncread(filename,'xco2_quality_flag');
pro_app=ncread(filename,'co2_profile_apriori');
akk=ncread(filename,'xco2_averaging_kernel');
%pl=ncread(filename,'pressure_levels'); make sure return the same name
%pl here. from D Baker, using psurf *sigma levels
uncc=ncread(filename,'xco2_uncertainty');
aprr=ncread(filename,'xco2_apriori');
co22=ncread(filename,'xco2');
pww=ncread(filename,'pressure_weight');
%load p_level.mat pl %calcuate beforehand and now load it here%

%GET PRESSURE LEVELS%
ps=ncread(filename,'psurf');
sl=ncread(filename,'sigma_levels');
pl=sl*ps'; clear sl ps; %hPa%




load ele.mat ele; %1 by 1 degree of elevation data, 1 for elevation <3000m; otherwise it's just 0%
%% this is the general word
for year =2015:2019,
    for month=1:12,
        for day=1:31,
            temp=zeros(1,1); %clear the value at each beginning of the day%
            k=0; %clear the value space%
            %%
            
            %%
            for i=1:length(datt(1,:)),
                
                %% used to filter the data with an altitude >4000m.
                
                %% 
                
                
                if datt(1,i)==year,
                    if datt(2,i)==month,
                        if datt(3,i)==day,
                            if latt(i)>-60, 
                                if latt(i)<60;%avoid the SNR issue%
                tedlon=round(lonn(i))+180; if tedlon==0; tedlon=1;end;
                tedlat=round(latt(i))+90;  if tedlat==0; tedlat=1;end;
                                if ele(tedlon,tedlat)==1,
                                %if rem(si(i),10)==1, %only take the land nadir%
                                 %if rem(si(i),10)==2, %only take the land glint%
                                 if (rem(si(i),10)==1)|(rem(si(i),10)==2), %nadir + glint land data% 
                                    k=k+1;
                                    temp(k,1)=i;
                                end
                            end
                        end
                    end
                end
                    end
                end
            end
            %% now we have the range to determine one day's data;
            
            %% below is to give the value for each day
            
            if length(temp)>5, %make sure we have at least more than five data points for the specific day%
                %define the variable here%
                length(temp)
                %%
                sid=zeros(1,1);
                dat=zeros(7,1);
                lon=zeros(1,1);
                lat=zeros(1,1);
                qff=zeros(1,1);
                profile_ap=zeros(20,1);
                ak=zeros(20,1);
                pls=zeros(20,1);
                unc=zeros(1,1);
                apr=zeros(1,1);
                co2=zeros(1,1);
                pw=zeros(20,1);
                %%
                
                %%we can also temp(end) to repalce the wordy code here;
                sid(1:length(temp))=si(temp);
                dat(:,1:length(temp))=datt(:,temp);
                lon(1:length(temp))=lonn(temp);
                lat(1:length(temp))=latt(temp);
                qff(1:length(temp))=qf(temp);
                profile_ap(:,1:length(temp))=pro_app(:,temp);
                unc(1:length(temp))=uncc(temp);
                apr(1:length(temp))=aprr(temp);
                co2(1:length(temp))=co22(temp);
                pw(:,1:length(temp))=pww(:,temp);
                pls(:,1:length(temp))=pl(:,temp);
                ak(:,1:length(temp))=akk(:,temp);
                %%
                
                
                %% BELOW is to read in
                %read in more than 1 variable%
                my_vardata = profile_ap; %co2_profile_apriori%
                my_vardata1 = dat;%date%
                my_vardata2 = sid;%sounding id%
                my_vardata3 = ak; %averaging kernel%
                my_vardata4 = pls;%pressure levels%
                my_vardata5 = lat;%latitude%
                my_vardata6 = lon; %longitude%
                my_vardata7 = unc;%xco2 uncertainty%
                my_vardata8 = apr;%xco2 apriopri%
                my_vardata9= co2; %xco2%
                my_vardata10 = pw;%pressure weight%
                
                % Open the netCDF file.
                ncid = netcdf.create('foo1.nc','NOCLOBBER');
                % Define the DIMENSIONS of the variable.
                dimid = netcdf.defDim(ncid,'sounding_id',length(sid));
                dimid1 = netcdf.defDim(ncid,'levels',20);
                dimid2 = netcdf.defDim(ncid,'epoch_dimension',7);
                
                % Define a new VARIALBE in the file.
                my_varID = netcdf.defVar(ncid,'co2_profile_apriori','double',[dimid1 dimid]);
                my_varID1 = netcdf.defVar(ncid,'date','double',[dimid2 dimid]);
                my_varID2 = netcdf.defVar(ncid,'sounding_id','double',dimid);
                my_varID3 = netcdf.defVar(ncid,'xco2_averaging_kernel','double',[dimid1 dimid]);
                my_varID4 = netcdf.defVar(ncid,'pressure_levels','double',[dimid1 dimid]);
                my_varID5 = netcdf.defVar(ncid,'latitude','double',dimid);
                my_varID6 = netcdf.defVar(ncid,'longitude','double',dimid);
                my_varID7 = netcdf.defVar(ncid,'xco2_uncertainty','double',dimid);
                my_varID8 = netcdf.defVar(ncid,'xco2_apriori','double',dimid);
                my_varID9 = netcdf.defVar(ncid,'xco2','double',dimid);
                my_varID10 = netcdf.defVar(ncid,'pressure_weight','double',[dimid1 dimid]);
                
                netcdf.endDef(ncid);
                % WRITE DATA to variable.
                netcdf.putVar(ncid,my_varID,my_vardata);
                netcdf.putVar(ncid,my_varID1,my_vardata1);
                netcdf.putVar(ncid,my_varID2,my_vardata2);
                netcdf.putVar(ncid,my_varID3,my_vardata3);
                netcdf.putVar(ncid,my_varID4,my_vardata4);
                netcdf.putVar(ncid,my_varID5,my_vardata5);
                netcdf.putVar(ncid,my_varID6,my_vardata6);
                netcdf.putVar(ncid,my_varID7,my_vardata7);
                netcdf.putVar(ncid,my_varID8,my_vardata8);
                netcdf.putVar(ncid,my_varID9,my_vardata9);
                netcdf.putVar(ncid,my_varID10,my_vardata10);
                
                % Verify that the variable was created.
                netcdf.close(ncid)
                
                
                
                if month<10,
                    month1=strcat('0',num2str(month));
                else
                    month1=num2str(month);
                end
                
                if day<10,
                    day1=strcat('0',num2str(day));
                else
                    day1=num2str(day);
                end
                
                
                
                movefile('foo1.nc',strcat('oco2_LtCO2_',num2str(year),month1,day1,'.nc'));
                %careful when using this command%
            end
            clearvars -except year month day si datt lonn latt qf pro_app akk pl uncc aprr co22 pww ele
        end
    end
end














