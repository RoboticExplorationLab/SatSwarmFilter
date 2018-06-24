function LunarData = DataReader(filename)
fileID = fopen(filename);
Data = textscan(fileID,'%*f %*c.%*c. %{yyyy-MMM-dd HH:mm:ss.SSSS}D %f %f %f %f %f %f','Delimiter',',','HeaderLines',23);
Time = Data{1};
X = Data{2}; %km
Y = Data{3};
Z = Data{4};
VX = Data{5}; %km/s
VY = Data{6};
VZ = Data{7};
Vec = [X,Y,Z,VX,VY,VZ];
LunarData = {Time, Vec};
fclose(fileID);