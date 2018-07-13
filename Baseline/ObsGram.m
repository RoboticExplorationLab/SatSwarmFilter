clear all
clc
beep on
N0 = 3; %Size of each target state
M = 2; %Number of total ships
N = 2*M*N0; %Size of total state
I = eye(N);
L = M*N0; %Size of observation

delT = minutes(10); %Timestep in minutes
lambda = 2; %UKF parameter
NOrbitSim = 1; %Number of orbits to be simulated (Only 1 for Obs Gram)
global muE muM Xtrue T_SIM LunarData
const = s3_constants;
muE = const.GM.EARTH/(1e3)^3; %km^3/s^2
muM = const.GM.MOON/(1e3)^3; %km^3/s^2

OE0 = [243729.554542e3,0.598335,51.975,61.110,247.345,0.000]; %Initial OE for hub in m
Xhub0 = s3_state_keptocart(OE0, 1)/1e3; %Initial ECI vector in km or km/s
a0 = (OE0(1)/1000); %Initial Semi Major Axis in km

Tapprox = seconds(2*pi*sqrt(a0^3/muE)); %Approximate obit period for estimating simulation time
% T0 = datetime(2024,01,25,22,30,25,520); %Initial Simulated Time
T0 = datetime(2024,01,25,22,30,25); %Initial Simulated Time
Tf = T0+Tapprox*NOrbitSim; %Final Simulated Time
T_SIM = (T0:delT:Tf).';

% Please download the lunar orbit data based on T0 and Tf from JPL
% ephemrides database: https://ssd.jpl.nasa.gov/horizons.cgi#top
DataName = 'LunaECI.txt'; %Data Input Name
LunarData = DataReader(DataName); %km & km/s
Q = ProcessNoise(M);
R = ObsNoise(M);

X0 = SampleGen(OE0,M); %Sample ECI vecotrs for all children but hub
X0{1} = Xhub0; %Input hub

Xtrue(:,1) = cell2mat(X0.');
mu0 = Xtrue(:,1) + sqrtm(Q)*randn(N,1);

Phi{1} = eye(N); %Initial state transition matrix
delTNum = seconds(delT);
for t = 1:numel(T_SIM)  %Compute true positions at each timestep
    T = T_SIM(t); %Extract real time
    W(:,t) = sqrtm(Q)*randn(N,1); %Process Noise
    
    Xtrue(:,t+1) = StateTrans(Xtrue(:,t),delT,T,LunarData)+W(:,t); %Compute true state for each timestep  
end
fprintf('Simulation Completed\n')
%% Compute Gramian
delT1 = duration(60*seconds(1));
T_Gram = (T0:delT1:Tf).'; 
delTNum = seconds(delT1);
k = 1;
for t = 1:numel(T_Gram)
    T = T_Gram(t);
    [At{t}, C{t}] = ACGenT(T);
    Phi{t+1} = exp(At{t}*delTNum)*Phi{t};
%     if At{t}
%         k = k+1
%     end
end
fprintf('State Transition Completed\n')
beep
% isSin = ~logical(det(W));
%%
%% Plotting
for i = 1:M
    plot3(Xtrue(3*(i-1)+1,:),Xtrue(3*(i-1)+2,:),Xtrue(3*(i-1)+3,:))
    hold on
end
plot3(LunarData{2}(:,1),LunarData{2}(:,2),LunarData{2}(:,3))