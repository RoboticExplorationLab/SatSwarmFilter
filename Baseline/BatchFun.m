function Diff = BatchFun(X0)
global T_SIM Y LunarData
delT = T_SIM(2) - T_SIM(1);
Xhat(:,1) = X0;
for t = 1:numel(T_SIM)  %Compute true positions at each timestep
    T = T_SIM(t); %Extract real time    
    Xhat(:,t+1) = StateTrans(Xhat(:,t),delT,T,LunarData); %Simulate each step
    Yhat(:,t) = Obs(Xhat(:,t+1)); %Compute observation from true state
end
Diff = Y - Yhat;
Diff = Diff(:);