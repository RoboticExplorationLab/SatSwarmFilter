clear
clc
N0 = 3; %Size of each target state
M = 4; %Number of total ships
N = 2*M*N0; %Size of total state
I = eye(N);
L = M*N0; %Size of observation


delT = 1;
lambda = 2;
Q = .1*I*delT;
R = .1;
mu0 = zeros(N,1);
sigma0 = .01*I;
Tf = 100;
Xtrue(:,1) = mu0;
tic
for t = 1:delT:Tf
    %Noise & Control
    W(:,t) = sqrtm(Q)*randn(N,1);
    V(:,t) = sqrtm(R)*randn(L,1); 
    u(:,t) = [1, sin(t)];
    %True State Simulation X(t+1) = f(X(t), u(t))+W(t)
    Xtrue(:,t+1) = StateTrans(Xtrue(:,t),delT)+W(:,t);
    Y(:,t) = Obs(Xtrue(:,t+1)) + V(:,t); %Compute observation from true state
    
    %Begin Prediction
    if t == 1
        Xsigma{t} = UT(mu0,sigma0);
    else
        Xsigma{t} = UT(mu_upd(:,t-1),sig_upd{t-1});
    end
    
    for j = 1:2*N+1
        Xpre_bar{t}(:,j) = StateTrans(Xsigma{t}(:,j), delT);
    end
    [mu_pre(:,t),sig_pre{t}] = invUT(Xpre_bar{t});
    sig_pre{t} = sig_pre{t} + Q;
    %End Prediction
    
    %Begin Update
    Xpre{t} = UT(mu_pre(:,t),sig_pre{t});
    for j = 1:2*N+1
        Ypre(:,j) = Obs(Xpre{t}(:,j));
    end     
    w0 = lambda/(lambda+N);
    wi = 1/(2*(lambda+N));
    
    Ypre_avg(:,t) = w0*Ypre(:,1);
    for j = 1:2*N
        Ypre_avg(:,t) = Ypre_avg(:,t) + wi*Ypre(:,j+1);         
    end
    sigyy{t} = w0*(Ypre(:,1)-Ypre_avg(:,t))*(Ypre(:,1)-Ypre_avg(:,t)).'+R;
    sigxy{t} = w0*(Xpre{t}(:,1)-mu_pre(:,t))*(Ypre(:,1)-Ypre_avg(:,t)).';
    for j = 1:2*N
        sigyy{t} = sigyy{t} + wi*(Ypre(:,j+1)-Ypre_avg(:,t))*(Ypre(:,j+1)-Ypre_avg(:,t)).';
        sigxy{t} = sigxy{t} + wi*(Xpre{t}(:,j+1)-mu_pre(:,t))*(Ypre(:,j+1)-Ypre_avg(:,t)).';
    end
    
    Kt = sigxy{t} * inv(sigyy{t});
    mu_upd(:,t) = mu_pre(:,t) + Kt * (Y(:,t) - Ypre_avg(:,t));
    sig_upd{t} = sig_pre{t} - Kt * sigxy{t}.';
    
    Xest(:,t) = sqrtm(sig_upd{t})*randn(N,1)+ mu_upd(:,t); %Estimation of X from update step
end
Xtrue(:,1) = []; %Remove the initial state
FinalTime = toc;
%%
StateLabel = ["p_{t}^{x}","p_{t}^{y}","\theta_{t}"];
TitleLabel = ["x Position with UKF","y Position with UKF","Heading Angle with UKF","Position with UKF"];
figure('Position', [100, 100, 1200, 750])
subplot(4,1,1)
plot(1:Tf,Xtrue(1,:),'-b',1:Tf,Xest(1,:),'--r','LineWidth',1)
title(TitleLabel(1),'FontSize',15)
xlabel('t','FontSize',15)
ylabel(StateLabel(1),'FontSize',15)
legend('True State','Estimation')

subplot(4,1,2)
plot(1:Tf,Xtrue(2,:),'-b',1:Tf,Xest(2,:),'--r','LineWidth',1)
title(TitleLabel(2),'FontSize',15)
xlabel('t','FontSize',15)
ylabel(StateLabel(2),'FontSize',15)
legend('True State','Estimation')

subplot(4,1,3)
plot(1:Tf,Xtrue(3,:),'-b',1:Tf,Xest(3,:),'--r','LineWidth',1)
title(TitleLabel(3),'FontSize',15)
xlabel('t','FontSize',15)
ylabel(StateLabel(3),'FontSize',15)
legend('True State','Estimation')

subplot(4,1,4)
plot(Xtrue(1,:),Xtrue(2,:),'-.b',Xest(1,:),Xest(2,:),'--*r','LineWidth',1)
title(TitleLabel(4),'FontSize',15)
xlabel(StateLabel(1),'FontSize',15)
ylabel(StateLabel(2),'FontSize',15)
legend('True State','Estimation')