clear; close; clc;
import matlab.io.hdfeos.*
import matlab.io.hdf4.*
%% Process TWV MYD05 Level 2 data%%
%%First lets get the lat lon data from Geolocation HDF file 
GEO_FILE_NAME='MYD03.A2018277.0740.061.2018277192140.hdf';
SWATH_NAME='MODIS_Swath_Type_GEO';
% Opening the HDF-EOS2 Swath File
file_id = sw.open(GEO_FILE_NAME, 'rdonly');
% Open swath
swath_id = sw.attach(file_id, SWATH_NAME);
% Read lat and lon dataset.
lon = sw.readField(swath_id, 'Longitude', [], [], []);
lat = sw.readField(swath_id, 'Latitude', [], [], []);
% Detach Swath object.
sw.detach(swath_id);
% Close file.
sw.close(file_id);
lon=double(lon);
lat=double(lat);
Lon=reshape(lon,[size(lon,1)*size(lon,2),1]);
Lat=reshape(lat,[size(lat,1)*size(lat,2),1]);
clear swath_id file_id SWATH_NAME;
%end
%%Data Field files%%
datafiles = dir(fullfile('/Users/pravash/Desktop/Paper_2/TWV/TCPW_NIR/', '*.hdf'));
FILE_NAME='MYD05_L2.A2018277.0740.061.2018277195948.hdf';
%filename = FILE_NAME.name;
newstr=split(FILE_NAME,".");
yeardate1=newstr(2,1);
yeardate=erase(yeardate1, 'A');
year=extractBefore(yeardate,5);
year=str2double(year);
date=extractAfter(yeardate,4);
date=str2double(date);
%date=datetime(date, 'convertfrom','juliandate','Format','yyy-MM-dd');
time=newstr(3,1);
hour=extractBefore(time,3);
min =extractAfter(time,2);
hour=str2double(hour);
min=str2double(min);
year=repmat(year,size(lat,1),size(lat,2));
year=reshape(year,[size(lat,1)*size(lat,2),1]);
date=repmat(date,size(lat,1),size(lat,2));
date=reshape(date,[size(lat,1)*size(lat,2),1]);
hour=repmat(hour,size(lat,1),size(lat,2));
hour=reshape(hour,[size(lat,1)*size(lat,2),1]);
min=repmat(min,size(lat,1),size(lat,2));
min=reshape(min,[size(lat,1)*size(lat,2),1]);
datelocation=[year,date,hour,min Lon,Lat];
clear date hour min newstr time year yeardate yeardate1; 
%% Now since date and location data are tabulated at 1kmX1km%%%
%% We can now extract the geophysical parameter and QA files for this granule%%
%% the parameter we are interested to extract are:
        %%Total Column Precipitable Water Vapor - Near Infrared Retrieval"
        %%Cloud_Mask_QA
        %%Quality_Assurance_Near_Infrared
file_id = sw.open(FILE_NAME, 'rdonly');
SWATH_NAME='mod05';
swath_id = sw.attach(file_id, SWATH_NAME);
DATAFIELD_NAME='Water_Vapor_Near_Infrared';
QA_ASSurance1='Cloud_Mask_QA'; %%QA1
QA_ASSurance2='Quality_Assurance_Near_Infrared'; %%QA2
%%%Reading the datafield%%%%
data1 = sw.readField(swath_id, DATAFIELD_NAME, [], [], []);
data=double(data1);
QA1 = sw.readField(swath_id, QA_ASSurance1, [], [], []);
QA1 = double(QA1);
QA2 = sw.readField(swath_id, QA_ASSurance2, [], [], []);
QA2=squeeze(QA2);
QA2=double(QA2);
%%Lets see some details of the dataset and then proceed 
%%1.adding offset and multiplying the scale factor%%
%%2. Fill value replaced
%%3.Work with QA1 and QA2
% Reading attributes from the data field
SD_id = sd.start(FILE_NAME, 'rdonly');
sds_index = sd.nameToIndex(SD_id, DATAFIELD_NAME);
sds_id = sd.select(SD_id, sds_index);
% get the long name of the data field
long_name_index = sd.findAttr(sds_id, 'long_name');
long_name = sd.readAttr(sds_id, long_name_index);
% Reading units from the data field
units_index = sd.findAttr(sds_id, 'unit');
units = sd.readAttr(sds_id, units_index);
% Reading filledValue from the data field
fillvalue_index = sd.findAttr(sds_id, '_FillValue');
fillvalue = sd.readAttr(sds_id, fillvalue_index);
% Reading scale_factor from the data field
scale_index = sd.findAttr(sds_id, 'scale_factor');
scale = sd.readAttr(sds_id, scale_index);
scale = double(scale);
% Reading add_offset from the data field
offset_index = sd.findAttr(sds_id, 'add_offset');
offset= sd.readAttr(sds_id, offset_index);
offset = double(offset);
% Reading valid_range from the data field
range_index = sd.findAttr(sds_id, 'valid_range');
range = sd.readAttr(sds_id, range_index);
% Terminate access to the corresponding data set
sd.endAccess(sds_id);
% Closing the File
sd.close(SD_id);
% Replacing the filled value with NaN
data(data==fillvalue) = NaN;
data(data > double(range(2))) = NaN;
data(data < double(range(1))) = NaN;
% Multiplying scale and adding offset, the equation is scale *(data-offset).
data = scale*(data-offset);
%%%Now lets get them all into one data frame%%%
data=reshape(data,[size(data1,1)*size(data1,2),1]);
QA1=reshape(QA1,[size(QA1,1)*size(QA1,2),1]);
QA2=reshape(QA2,[size(QA2,1)*size(QA2,2),1]);
TWV_NOQC =[datelocation data QA1 QA2];
clearvars -except TWV_NOQC
%% We will first figure out  Cloud_Mask_QA
%% let us convert the QC Cloud Mask from decimal to binary
QC1  = string(dec2bin(TWV_NOQC(:,8), 8));
%% CM_QA is a 1 byte dat with 8 bits where bits are read 
%% the bits will always be numbered from right (bit index #1) to left (bit index #8)
%% in the combination of (8,(7,6),5,4,3,(2,1))
%% seperating bit by bit as per (Cloud Mask QA)
CMQA(:,1)=extract(QC1(:,1),8);
CMQA(:,2)=extractBetween(QC1(:,1),6,7);
CMQA(:,3)=extract(QC1(:,1),5);
CMQA(:,4)=extract(QC1(:,1),4);
CMQA(:,5)=extract(QC1(:,1),3);
CMQA(:,6)=extractBetween(QC1(:,1),1,2);
%% Converting each column from binary to decimal and referring lookup table
CMQA1(:,1)=bin2dec(CMQA(:,1)); %% CM_flagstatus (1)
CMQA1(:,2)=bin2dec(CMQA(:,2)); %% CM_cloudiness (3)
CMQA1(:,3)=bin2dec(CMQA(:,3)); %% Night_day (1)
CMQA1(:,4)=bin2dec(CMQA(:,4)); %%Sunlight_flag (1)
CMQA1(:,5)=bin2dec(CMQA(:,4)); %%Ice flag (1)
CMQA1(:,6)=bin2dec(CMQA(:,6)); %%Surface Type Flag (Fill for ocean) (0)
%% Next, will figure out Quality_Assurance_Near_Infrared(QA_NIR)
QC2=string(dec2bin(TWV_NOQC(:,9), 8));
%% QA_NIR is a 1 byte dat with 8 bits where bits are read 
%% the bits will always be numbered from right (bit index #1) to left (bit index #8)
%% in the combination of (8,(7,6,5),(4,3),(2,1))
%% seperating bit by bit as per (QA_NIR)
NIRQA(:,1)=extract(QC2(:,1),8);
NIRQA(:,2)=extractBetween(QC2(:,1),5,7);
NIRQA(:,3)=extractBetween(QC2(:,1),3,4);
NIRQA(:,4)=extractBetween(QC2(:,1),1,2);
%% Converting each column from binary to decimal and referring lookup table
NIRQA1(:,1)=bin2dec(NIRQA(:,1)); %% Total Precpitable water Usefulness (1)
NIRQA1(:,2)=bin2dec(NIRQA(:,2)); %% Total Precpitable water Confidence (2 and 3)
NIRQA1(:,3)=bin2dec(NIRQA(:,3)); %% Inversion method (0 and 1)
NIRQA1(:,4)=bin2dec(NIRQA(:,4)); %%Surface type (0)
%% Combining all the dataset and processed Quality Assurance values
TWV_NOQC1=[TWV_NOQC CMQA1 NIRQA1];
%% Subset only lat-lon ranges that are within within and around Bangladesh boundar %%
idx1 =  TWV_NOQC1(:,5) >=87.4 & TWV_NOQC1(:,5) <=93.0  &  TWV_NOQC1(:,6) >= 20.75 & TWV_NOQC1(:,6)<= 27.0;
TWV_NOQC1=TWV_NOQC1(idx1, :);
%% Incorporating all the quality flags 
idx2 = TWV_NOQC1(:,10)>0 & TWV_NOQC1(:,11)>2 & TWV_NOQC1(:,12)>0 & TWV_NOQC1(:,13)>0 &  TWV_NOQC1(:,14)>0 & TWV_NOQC1(:,15)>0; ...
    TWV_NOQC1(:,16)>0 & TWV_NOQC1(:,17)>1 & TWV_NOQC1(:,18)<2 & TWV_NOQC1(:,19)<1;
TCPW=TWV_NOQC1(idx2, :);
TCPW(:,8:19) = [];
%% look up table and QA GPlan for C61: https://modis-images.gsfc.nasa.gov/_docs/QA_Plan_C61_Master_2017_03_15.pdf 
