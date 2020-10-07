%% cc Milo Group BGU 2020 
% modified from Milo et al. Nature, 2014, 507, 210. 
%%

clc
clear all
close all

%% import data

% uncomment below for 2-sub rates
%[data,textdata,raw] = xlsread('datasetsLR','rate','A1:T12'); 

% uncomment below for catalyst set
[data,textdata,raw] = xlsread('datasetsLR','cat','A1:AA16'); 

%% define parameter set to be tested 

% uncomment below for 2-sub rates
%param = data(:,[1:(end-2)]);

% uncomment below for catalyst set
param = data(:,[1:(end-4)]);

deut = data(:,(end)); 

% for correlation matrix

% uncomment below for 2-sub rates
%corrParam = data(:,1:(end-2));
%nameParam = textdata(2,2:(end-2)); 

% uncomment below for catalyst set
corrParam = data(:,1:(end-4)); 
nameParam = textdata(2,2:(end-4));

%% normalize parameters and output by subtracting mean and deviding by standard deviation

stdtCat = nanstd(param);
meanCat = nanmean(param);
srCat1 = param-repmat(meanCat,size(param,1),1);
srParam = srCat1./repmat(stdtCat,size(param,1),1);

stdtOut = nanstd(deut);
meanOut = nanmean(deut);
srOut1 = deut-repmat(meanOut,size(deut,1),1);
srOut = srOut1./repmat(stdtOut,size(deut,1),1);

stdtCorrCat = nanstd(corrParam); % for correlation matrix
meanCorrCat = nanmean(corrParam); % for correlation matrix
srCorrCat1 = corrParam-repmat(meanCorrCat,size(corrParam,1),1); % for correlation matrix
srCorrParam = srCorrCat1./repmat(stdtCorrCat,size(corrParam,1),1); % for correlation matrix

%% define X and Y

Y = deut; 

% uncomment for model with normalized parameters 
X = srParam(:,:);

% uncomment for model with parameters not normalized
%X = param(:,:);

%% parameter correlation matrix 
% uses Correlations - modified from - Guo et al. ACS Catal. 2017, 7, 4144.
% uses brewermap - % (c) 2015 Stephen Cobeldick & (c) 2002 Cynthia Brewer, 
% Mark Harrower, and The Pennsylvania State University.

Correlations(srCorrParam, nameParam);

%% model development (Jane Fonda: Fit & Strong Level 1)

% 4 automated modeling algorithms: 

% add terms according to a p-value threshold (the default entrance 
% tolerance is 0.05) then remove terms from a linear fit according to a 
% p-value threshold (the default exit tolerance of 0.1) & iterate until 
% statistically significant model

FitForward = LinearModel.stepwise(X, Y)

% remove terms starting from a linear fit (a constant and linear terms: 
% x1 + x2 + x3 etc.) according to a p-value threshold then add terms 
% according to a p-value threshold & iterate

FitBackL = LinearModel.stepwise(X, Y, 'linear')

% remove terms then add terms according to p-value starting from an 
% interaction model (a constant, linear terms & interaction terms: 
% x1 + x2 + x3 + x1:x2 + x1:x3 + x2:x3) & iterate

%FitBackI = LinearModel.stepwise(X, Y, 'interactions')

% remove terms then add terms according to p-value starting from a 
% quadratic model (a constant, linear terms & squared terms: 
% x1 + x2 + x3 + x1^2 + x2^2 + x3^2) & iterate

%FitQuad = LinearModel.stepwise(X, Y, 'purequadratic')

%% interactive modeling (Jane Fonda: Fit & Strong Level 2)
 
% stepwise regression (Step) from and linear regression (Fit) to a chosen 
% model. Iterate by replacing the parameters (pink) 

% uncomment below for 2-sub rates
%Step = LinearModel.stepwise(X, Y, 'y ~ x1 + x5')
%Fit = LinearModel.fit(X, Y, 'y ~ x1 + x5')

% uncomment below for catalyst set
Step = stepwiselm(X, Y, 'y ~ x14 + x22')
Fit = fitlm(X, Y, 'y ~ x14 + x22')

%% predict and plot results (measured vs predicted) 

figure;
Ypred = predict(Fit,X);
Ysy = [Y,Ypred];
axis equal;
plotregression(Ysy(1:end,1),Ysy(1:end,2),'Regression')

%% K-fold cross validation - uses K_Fold_CV_avgN - Guo et al. ACS Catal. 2017, 7, 4144.

% uncomment below for 2-sub rates
%[Q2,Qsq] = K_Fold_CV_avgN(X,Y,'y ~ 1 + x1 + x5',3,500);

% uncomment below for catalyst set
[Q2,Qsq] = K_Fold_CV_avgN(X,Y,'y ~ 1 + x14 + x22',3,500);

%% caluculate and plot leave one out

Xp = X;
Yp = Y;

for i = 1:size(Xp,1)

    valY = Yp(i);
    valX = Xp(i,:);

    if i == 1
        X = Xp(2:end,:);
        Y = Yp(2:end);
    elseif i == size(Xp,1)
        X = Xp(1:(size(Xp,1)-1),:);
        Y = Yp(1:(size(Yp,1)-1));
    else
        Y = [Yp(1:i-1);Yp(i+1:size(Yp,1))];
        X = [Xp(1:i-1,:);Xp(i+1:size(Xp,1),:)];
    end
    
% uncomment below for 2-sub rates
%    Fit = LinearModel.fit(X, Y, 'y ~ x1 + x5')
    
% uncomment below for catalyst set
    Fit = LinearModel.fit(X, Y, 'y ~ x14 + x22')

    Ypredloo = predict(Fit,valX);
    y(i).yval = valY;
    y(i).ypred = Ypredloo;
end

for i = 1:size(y,2)
    
    Ys(i,1) = y(i).yval;
    Ys(i,2) = y(i).ypred; 
end
 
X = Xp;
Y = Yp;

figure;
axis equal;
plotregression(Ys(1:end,1),Ys(1:end,2),'Regression')




