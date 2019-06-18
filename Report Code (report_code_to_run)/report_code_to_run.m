%%% RUN THIS FILE INORDER TO GET SIMULATIONS FROM REPORT
%%% PLEASE NOTE THAT THE TABLES WOULD BE REQUIRED TO RUN SEPARATELY

%%% The starter section below is from: 
%%% http://www.columbia.edu/~fy2158/software/ROAD.html
p = 1000;     % number of variables
n = 300;      % number of observations
s0 = 10;      % number of nonzero mean differences
rho = 0.5;    % pairwise correlation coefficient
randSeed = 1; % setup the random seed to let the program repeatable
mu1 = zeros(p,1);
mu2 = zeros(p,1);
mu2(1:s0) = 1;
Sigma = eqcor(p,rho);
rand('state', randSeed);
randn('state',randSeed);
nTrain  = n;
nTest   = n;
Y1Train = mvnrnd(repmat(mu1',nTrain,1),Sigma);
Y2Train = mvnrnd(repmat(mu2',nTrain,1),Sigma);
Y1Test  = mvnrnd(repmat(mu1',nTest,1),Sigma);
Y2Test  = mvnrnd(repmat(mu2',nTest,1),Sigma);
x = [Y1Train;Y2Train];
y = [zeros(n,1);ones(n,1)];
xtest = [Y1Test;Y2Test];
ytest = [zeros(n,1);ones(n,1)];
[ROADfit] = roadBatch(x, y, xtest, ytest, 0, 0); % ROAD
[sROAD1fit] = roadBatch(x, y, xtest, ytest, 0, 1); % Screening-based ROAD version 1 (S-ROAD1)
[sROAD2fit] = roadBatch(x, y, xtest, ytest, 0, 2); % Screening-based ROAD version 2 (S-ROAD2)

%%% -----------------------------------------------------------------------
%%% Only the code below this line is our own code
%%% -----------------------------------------------------------------------
[DROADfit] = roadBatch(x, y, xtest, ytest, 1, 0); % DROAD

% Plotting out the different solution paths for ROAD (Simulated Data)
figure(1)
plot((1:100), ROADfit.wPath')
title("ROAD Path")
xlabel("Penalty Parameter Index")
ylabel("Coefficient")

% Plotting out the different solution paths for DROAD, SROAD1, SROAD2
% (Simulated Data)
figure(2)
plot((1:100), DROADfit.wPath')
title("DROAD Path")
xlabel("Penalty Parameter Index")
ylabel("Coefficient")

figure(3)
plot((1:100), sROAD1fit.wPath')
title("sRoad1 Path")
xlabel("Penalty Parameter Index")
ylabel("Coefficient")


figure(4)
plot((1:100), sROAD2fit.wPath')
title("sRoad2 Path")
xlabel("Penalty Parameter Index")
ylabel("Coefficient")


% Simulating the solution paths of the real data
 
% Data from: https://github.com/ramhiser/datamicroarray/blob/master/data/golub.RData
train = csvread('GOLUBtrain.csv',1,1);
test = csvread('GOLUBtest.csv', 1,1);

x_train = train(:,1:7129);
y_train = train(:,7130);

x_test = test(:,1:7129);
y_test = test(:,7130);

[gROADfit] = roadBatch(x_train, y_train, x_test, y_test, 0, 0); % ROAD

[gsROAD1fit] = roadBatch(x_train, y_train, x_test, y_test, 0, 1); % (S-ROAD1)

[gsROAD2fit] = roadBatch(x_train, y_train, x_test, y_test, 0, 2); % (S-ROAD2)


figure(5)
plot((1:100), gROADfit.wPath')
title("Road Path on Real Data")
xlabel("Penalty Parameter Index")
ylabel("Coefficient")

figure(6)
plot((1:100), gsROAD1fit.wPath')
title("sRoad1 Path on Real Data")
xlabel("Penalty Parameter Index")
ylabel("Coefficient")

figure(7)
plot((1:100), gsROAD2fit.wPath')
title("sRoad2 Path on Real Data")
xlabel("Penalty Parameter Index")
ylabel("Coefficient")

% Classification Error for the Leukaemia Data
Method = {'ROAD';'S-ROAD1';'S-ROAD2'};
Training_Error = [gROADfit.trainError;gsROAD1fit.trainError;gsROAD2fit.trainError];
Testing_Error = [gROADfit.testError;gsROAD1fit.testError;gsROAD2fit.testError];
rs0 = gROADfit.lamind;
rs1 = gsROAD1fit.lamind;
rs2 = gsROAD2fit.lamind;
Number_of_Selected_Genes = [nnz(gROADfit.wPath(:,rs0));nnz(gsROAD1fit.wPath(:,rs1));nnz(gsROAD2fit.wPath(:,rs2))];
% Table used in report
T = table(Method,Training_Error,Testing_Error,Number_of_Selected_Genes);
T;

% Displaying how Correlation Affects the Number of Non-zero coefficients in
% ROAD methods

ROAD = zeros(10,1);
SROAD1 = zeros(10,1);
SROAD2 = zeros(10,1);
DROAD = zeros(10,1);
rho_val = zeros(10,1);
i = 0;
for new_rho = 0:0.1:0.9
    i = i + 1;
    Sigma = eqcor(p,new_rho);
    rand('state', randSeed);
    randn('state',randSeed);
    nTrain  = n;
    nTest   = n;
    Y1Train = mvnrnd(repmat(mu1',nTrain,1),Sigma);
    Y2Train = mvnrnd(repmat(mu2',nTrain,1),Sigma);
    Y1Test  = mvnrnd(repmat(mu1',nTest,1),Sigma);
    Y2Test  = mvnrnd(repmat(mu2',nTest,1),Sigma);
    x = [Y1Train;Y2Train];
    y = [zeros(n,1);ones(n,1)];
    xtest = [Y1Test;Y2Test];
    ytest = [zeros(n,1);ones(n,1)];

    [cROADfit] = roadBatch(x, y, xtest, ytest, 0, 0);
    [csROAD1fit] = roadBatch(x, y, xtest, ytest, 0, 1);
    [csROAD2fit] = roadBatch(x, y, xtest, ytest, 0, 2);
    [cDROADfit] = roadBatch(x, y, xtest, ytest, 1, 0);
    
    c0 = cROADfit.lamind;   
    c1 = csROAD1fit.lamind;
    c2 = csROAD2fit.lamind;
    c3 = cDROADfit.lamind;
    rho_val(i) = new_rho;
    ROAD(i) = nnz(cROADfit.wPath(:,c0));
    SROAD1(i) = nnz(csROAD1fit.wPath(:,c1));
    SROAD2(i) = nnz(csROAD2fit.wPath(:,c2));
    DROAD(i) = nnz(cDROADfit.wPath(:,c3));
end

% Table used in Report
methods_Table = table(rho_val, ROAD, SROAD1, SROAD2, DROAD);
methods_Table;