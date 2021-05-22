%% Laboratory 5: Hurst exponent
%% Jairo Viola
clear all; close all;clc;
%% part 1: Generating fGn with known Hurst parameters
%%
% 
%  The fractional gaussian noise function was changed to this one
%  https://www.mathworks.com/matlabcentral/mlc-downloads/downloads/submissions/29686/versions/1/previews/ffGn.m/index.html
%  due to the result of the other one is not consistent whe the hurst
%  parameter is calculated by any of the proposed estimators in the
%  toolbox.
%  
%% Generating FGN 
hurst=0.1:0.1:0.9
for i=1:length(hurst)
    hurst1=ffGn(4096,hurst(i),1,0);
    subplot(3,3,i)
    plot(hurst1)
    titleText=strcat("FGN for H=",string(hurst(i)));
    title(titleText)
    set(gca,'FontSize', 10);
    xlim([0 4000])
end




%% theoretical and estimated hurst exponents on nominal conditions
%hurst parameter
H=0.05:0.05:1;
%fGN generation
sigma=1;  %variance
n=4096;   %number of samples  

%number of timeseries generated
N=1:round(100/length(H)):100;    
M=n/2;
force=1;

%Hurst parameter plot calculation
for i=1:length(H)
%     M=20;
    for k=1:N
        aux1=ffGn(1000,H(i),1,0);
        H1(k)=hurst_estimate(aux1','RS');
        H2(k)=hurst_estimate(aux1','aggvar');
        H3(k)=hurst_estimate(aux1','boxper');
        H4(k)=hurst_estimate(aux1','diffvar');
        H5(k)=hurst_estimate(aux1','absval');
        H6(k)=hurst_estimate(aux1','peng');
    end
    H1Hist(i)=mean(H1); H2Hist(i)=mean(H2); H3Hist(i)=mean(H3);
    H4Hist(i)=mean(H4); H5Hist(i)=mean(H5); H6Hist(i)=mean(H6);
    
end

Horig=H';
H_RS=H1Hist'; H_aggvar=H2Hist'; H_boxper=H3Hist';
H_diffvar=H4Hist'; H_absval=H5Hist';    H_peng=H6Hist';
T = table(Horig, H_RS, H_aggvar, H_boxper, H_diffvar, H_absval,H_peng);

%estimation error calculation
TNew=table2array(T);
TNew(isnan(TNew))=0;
[row col]=size(TNew);
TError=TNew;
for j=2:col
   TError(:,j)=abs(((TNew(:,j)-TNew(:,1)))./TNew(:,1))*100;
end
TError
figure()
for i=2:col
    plot(TError(:,1),TError(:,i))
    legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
    hold on
end

figure()
subplot(2,1,1)
bar(TError(1:8),TError(1:8,:))
legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
title('Estimation error for 0.05<H<0.4')
xlabel('Hurst exponent')
ylabel('error')
set(gca,'FontSize', 12);

subplot(2,1,2)
bar(TError(9:20),TError(9:end,:))
legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
title('Estimation error for 0.45<H<1')
xlabel('Hurst exponent')
ylabel('error')
set(gca,'FontSize', 12);

%%
% 
%  On the next figure it is posible to see the estimation of hurst
%  parameters employing six different estimation criteria:
%  RS, aggregated value, periodogram, difference variance method, 
%  absuolute moment, and regression residuals method.
%  
% It can be observed that each hurst parameter estimator has itown
% advantages at different LRD levels. Therefore the best option to evaluate
% the LDR of a time series is employing multiples estimators to decide
% which is the best estimator for the data.
%

figure()
plot(N,H)
hold on
loglog(N,(H1Hist),'-*')
loglog(N,(H2Hist),'-*')
loglog(N,(H3Hist),'-*')
loglog(N,(H4Hist),'-*')
loglog(N,(H5Hist),'-*')
loglog(N,(H6Hist),'-*')
% loglog(N,(H7Hist))
legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
xlabel('Number of Hurst exponents calculated')
ylabel('Hurst exponent')
title('Hurst exponent for multiple Estimators')
set(gca,'FontSize', 12);
%% theoretical and estimated hurst exponents with 30dB gaussian noise
clear all;close all;clc
%hurst parameter
H=0.05:0.05:1;
%fGN generation with variance 1 and 500 samples
sigma=1;
n=4096;
%number of timeseries generated per hurst exponent
N=1:round(100/length(H)):100;    %number of timeseries generated
M=n/2;
force=1;
%Hurst parameter plot calculation
for i=1:length(H)
%     M=20;
    for k=1:length(N)
        aux1=ffGn(1000,H(i),1,0);
        aux1 = awgn(aux1,30);
     
        H1(k)=hurst_estimate(aux1','RS');
        H2(k)=hurst_estimate(aux1','aggvar');
        H3(k)=hurst_estimate(aux1','boxper');
        H4(k)=hurst_estimate(aux1','diffvar');
        H5(k)=hurst_estimate(aux1','absval');
        H6(k)=hurst_estimate(aux1','peng');
    end
    H1Hist(i)=mean(H1);
    H2Hist(i)=mean(H2);
    H3Hist(i)=mean(H3);
    H4Hist(i)=mean(H4);
    H5Hist(i)=mean(H5);
    H6Hist(i)=mean(H6);
    
end
Horig=H';
H_RS=H1Hist'; H_aggvar=H2Hist'; H_boxper=H3Hist';
H_diffvar=H4Hist'; H_absval=H5Hist';    H_peng=H6Hist';
T = table(Horig, H_RS, H_aggvar, H_boxper, H_diffvar, H_absval,H_peng)

Horig=H';
H_RS=H1Hist'; H_aggvar=H2Hist'; H_boxper=H3Hist';
H_diffvar=H4Hist'; H_absval=H5Hist';    H_peng=H6Hist';
T = table(Horig, H_RS, H_aggvar, H_boxper, H_diffvar, H_absval,H_peng);

%estimation error calculation
TNew=table2array(T);
TNew(isnan(TNew))=0;
[row col]=size(TNew);
TError=TNew;
for j=2:col
   TError(:,j)=abs(((TNew(:,j)-TNew(:,1)))./TNew(:,1))*100;
end
TError
figure()
for i=2:col
    plot(TError(:,1),TError(:,i))
    legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
    hold on
end

figure()
subplot(2,1,1)
bar(TError(1:8),TError(1:8,:))
legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
title('Estimation error for 0.05<H<0.4')
xlabel('Hurst exponent')
ylabel('error')
set(gca,'FontSize', 12);

subplot(2,1,2)
bar(TError(9:20),TError(9:end,:))
legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
title('Estimation error for 0.45<H<1')
xlabel('Hurst exponent')
ylabel('error')
set(gca,'FontSize', 12);


%%
% 
% The following figure shows the robustness of the Hurst parameter
% estimators in the presence of white gaussian noise with 30 dB of
% intensity. It can be observed that the regression residuals method is
% more robust on the estimation of the hurst parameter in the presence of
% random noise but fails as H=1. The same behavior repeats with the other
% estimators evaluated.
%

figure()
plot(N,H)
hold on
loglog(N,(H1Hist),'-*')
loglog(N,(H2Hist),'-*')
loglog(N,(H3Hist),'-*')
loglog(N,(H4Hist),'-*')
loglog(N,(H5Hist),'-*')
loglog(N,(H6Hist),'-*')
% loglog(N,(H7Hist))
legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
title('Hurst exponent with 30db SNR white gaussian noise')
legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
xlabel('Number of Hurst exponents calculated')
ylabel('Hurst exponent')
set(gca,'FontSize', 12);
%% theoretical and estimated hurst exponents with 30dB alpha stable noise
clear all;close all;clc
%hurst parameter
H=0.05:0.05:1;
%fGN generation with variance 1 and 500 samples
sigma=1;
n=4096;
%number of timeseries generated per hurst exponent
N=1:round(100/length(H)):100;    %number of timeseries generated
M=n/2;
force=1;
%Hurst parameter plot calculation
for i=1:length(H)
%     M=20;
    for k=1:length(N)
        aux1=ffGn(1000,H(i),1,0);
        aux1 = aux1+stblrnd(0.9,0,1,0.1);
     
        H1(k)=hurst_estimate(aux1','RS');
        H2(k)=hurst_estimate(aux1','aggvar');
        H3(k)=hurst_estimate(aux1','boxper');
        H4(k)=hurst_estimate(aux1','diffvar');
        H5(k)=hurst_estimate(aux1','absval');
        H6(k)=hurst_estimate(aux1','peng');
    end
    H1Hist(i)=mean(H1);
    H2Hist(i)=mean(H2);
    H3Hist(i)=mean(H3);
    H4Hist(i)=mean(H4);
    H5Hist(i)=mean(H5);
    H6Hist(i)=mean(H6);
    
end
Horig=H';
H_RS=H1Hist'; H_aggvar=H2Hist'; H_boxper=H3Hist';
H_diffvar=H4Hist'; H_absval=H5Hist';    H_peng=H6Hist';
T = table(Horig, H_RS, H_aggvar, H_boxper, H_diffvar, H_absval,H_peng)

Horig=H';
H_RS=H1Hist'; H_aggvar=H2Hist'; H_boxper=H3Hist';
H_diffvar=H4Hist'; H_absval=H5Hist';    H_peng=H6Hist';
T = table(Horig, H_RS, H_aggvar, H_boxper, H_diffvar, H_absval,H_peng);

%estimation error calculation
TNew=table2array(T);
TNew(isnan(TNew))=0;
[row col]=size(TNew);
TError=TNew;
for j=2:col
   TError(:,j)=abs(((TNew(:,j)-TNew(:,1)))./TNew(:,1))*100;
end
TError
figure()
for i=2:col
    plot(TError(:,1),TError(:,i))
    legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
    hold on
end

figure()
subplot(2,1,1)
bar(TError(1:8),TError(1:8,:))
legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
title('Estimation error for 0.05<H<0.4')
xlabel('Hurst exponent')
ylabel('error')
set(gca,'FontSize', 12);

subplot(2,1,2)
bar(TError(9:20),TError(9:end,:))
legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
title('Estimation error for 0.45<H<1')
xlabel('Hurst exponent')
ylabel('error')
set(gca,'FontSize', 12);

%%
% 
% The following figure shows the robustness of the Hurst parameter
% estimators in the presence of symetric alpha stable random noise with 30 dB of
% intensity and $alpha=0.9$. In this case the presence of alpha stable noise shows 
% that still the residual regression method is capable to obtain a good
% response in the presence of random nois. Also THe absolute value and
% aggregated value Hurst parameters exhibit a good performance for negative
% LRD timeseries.
%
%
figure()
plot(N,H)
hold on
loglog(N,(H1Hist),'-*')
loglog(N,(H2Hist),'-*')
loglog(N,(H3Hist),'-*')
loglog(N,(H4Hist),'-*')
loglog(N,(H5Hist),'-*')
loglog(N,(H6Hist),'-*')
% loglog(N,(H7Hist))
legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
title('Hurst exponent with 30db SNR \alpha stable noise')
legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
xlabel('Number of Hurst exponents calculated')
ylabel('Hurst exponent')
set(gca,'FontSize', 12);


%% part 2: Decimation effects over H parameter estimation
%%
% 
% In this section a subsampling process is performed to the fGn timeseries
% in order to analyze its effects on the hurst parameter estimation
%
%  
%% theoretical and estimated hurst exponents on nominal conditions
%hurst parameter
H=0.05:0.05:1;
%fGN generation with variance 1 and 500 samples
sigma=1;
n=4096;
%number of timeseries generated per hurst exponent
N=1:round(100/length(H)):100;    %number of timeseries generated
M=n/2;
force=1;
%Hurst parameter plot calculation
for i=1:length(H)
%     M=20;
    for k=1:length(N)
        auxO=ffGn(n,H(i),1,0);
        aux1 = decimate(auxO,2);
        H1(k)=hurst_estimate(aux1','RS');
        H2(k)=hurst_estimate(aux1','aggvar');
        H3(k)=hurst_estimate(aux1','boxper');
        H4(k)=hurst_estimate(aux1','diffvar');
        H5(k)=hurst_estimate(aux1','absval');
        H6(k)=hurst_estimate(aux1','peng');
    end
    H1Hist(i)=mean(H1);
    H2Hist(i)=mean(H2);
    H3Hist(i)=mean(H3);
    H4Hist(i)=mean(H4);
    H5Hist(i)=mean(H5);
    H6Hist(i)=mean(H6);
    
end

Horig=H';
H_RS=H1Hist'; H_aggvar=H2Hist'; H_boxper=H3Hist';
H_diffvar=H4Hist'; H_absval=H5Hist';    H_peng=H6Hist';
T = table(Horig, H_RS, H_aggvar, H_boxper, H_diffvar, H_absval,H_peng)

Horig=H';
H_RS=H1Hist'; H_aggvar=H2Hist'; H_boxper=H3Hist';
H_diffvar=H4Hist'; H_absval=H5Hist';    H_peng=H6Hist';
T = table(Horig, H_RS, H_aggvar, H_boxper, H_diffvar, H_absval,H_peng);

%estimation error calculation
TNew=table2array(T);
TNew(isnan(TNew))=0;
[row col]=size(TNew);
TError=TNew;
for j=2:col
   TError(:,j)=abs(((TNew(:,j)-TNew(:,1)))./TNew(:,1))*100;
end
TError
figure()
for i=2:col
    plot(TError(:,1),TError(:,i))
    legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
    hold on
end

figure()
subplot(2,1,1)
bar(TError(1:8),TError(1:8,:))
legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
title('Estimation error for 0.05<H<0.4')
xlabel('Hurst exponent')
ylabel('error')
set(gca,'FontSize', 12);

subplot(2,1,2)
bar(TError(9:20),TError(9:end,:))
legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
title('Estimation error for 0.45<H<1')
xlabel('Hurst exponent')
ylabel('error')
set(gca,'FontSize', 12);




%% 
%
% From the table and figure  it can be observed that the 
% peng and absolute value estimators obtaines the best estimation
% for the Hurst parameter.
%

Horig=H';
H_RS=H1Hist'; H_aggvar=H2Hist'; H_boxper=H3Hist';
H_diffvar=H4Hist'; H_absval=H5Hist';    H_peng=H6Hist';
T = table(Horig, H_RS, H_aggvar, H_boxper, H_diffvar, H_absval,H_peng)

figure()
plot(N,H)
hold on
loglog(N,(H1Hist),'-*')
loglog(N,(H2Hist),'-*')
loglog(N,(H3Hist),'-*')
loglog(N,(H4Hist),'-*')
loglog(N,(H5Hist),'-*')
loglog(N,(H6Hist),'-*')
% loglog(N,(H7Hist))
legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
title('Hurst exponent for multiple Estimators')
xlabel('Number of Hurst exponents calculated decimated')
ylabel('Hurst exponent')
set(gca,'FontSize', 12);



%% theoretical and estimated hurst exponents with 30dB gaussian noise
clear all;close all;clc
%hurst parameter
H=0.05:0.05:1;
%fGN generation with variance 1 and 500 samples
sigma=1;
n=4096;
%number of timeseries generated per hurst exponent
N=1:round(100/length(H)):100;    %number of timeseries generated
M=n/2;
force=1;
%Hurst parameter plot calculation
for i=1:length(H)
%     M=20;
    for k=1:length(N)
        auxO=ffGn(n,H(i),1,0);
        aux1 = decimate(auxO,2);
        aux1 = awgn(aux1,30);
     
        H1(k)=hurst_estimate(aux1','RS');
        H2(k)=hurst_estimate(aux1','aggvar');
        H3(k)=hurst_estimate(aux1','boxper');
        H4(k)=hurst_estimate(aux1','diffvar');
        H5(k)=hurst_estimate(aux1','absval');
        H6(k)=hurst_estimate(aux1','peng');
    end
    H1Hist(i)=mean(H1);
    H2Hist(i)=mean(H2);
    H3Hist(i)=mean(H3);
    H4Hist(i)=mean(H4);
    H5Hist(i)=mean(H5);
    H6Hist(i)=mean(H6);
    
end

Horig=H';
H_RS=H1Hist'; H_aggvar=H2Hist'; H_boxper=H3Hist';
H_diffvar=H4Hist'; H_absval=H5Hist';    H_peng=H6Hist';
T = table(Horig, H_RS, H_aggvar, H_boxper, H_diffvar, H_absval,H_peng)

Horig=H';
H_RS=H1Hist'; H_aggvar=H2Hist'; H_boxper=H3Hist';
H_diffvar=H4Hist'; H_absval=H5Hist';    H_peng=H6Hist';
T = table(Horig, H_RS, H_aggvar, H_boxper, H_diffvar, H_absval,H_peng);

%estimation error calculation
TNew=table2array(T);
TNew(isnan(TNew))=0;
[row col]=size(TNew);
TError=TNew;
for j=2:col
   TError(:,j)=abs(((TNew(:,j)-TNew(:,1)))./TNew(:,1))*100;
end
TError
figure()
for i=2:col
    plot(TError(:,1),TError(:,i))
    legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
    hold on
end

figure()
subplot(2,1,1)
bar(TError(1:8),TError(1:8,:))
legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
title('Estimation error for 0.05<H<0.4')
xlabel('Hurst exponent')
ylabel('error')
set(gca,'FontSize', 12);

subplot(2,1,2)
bar(TError(9:20),TError(9:end,:))
legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
title('Estimation error for 0.45<H<1')
xlabel('Hurst exponent')
ylabel('error')
set(gca,'FontSize', 12);


%%
%
%In the presence of random noise with 30dB SNR, the 
%residual covariance and the diff var methods returns the best hurst
%estimatior performance.
%

Horig=H';
H_RS=H1Hist'; H_aggvar=H2Hist'; H_boxper=H3Hist';
H_diffvar=H4Hist'; H_absval=H5Hist';    H_peng=H6Hist';
T = table(Horig, H_RS, H_aggvar, H_boxper, H_diffvar, H_absval,H_peng)


figure()
plot(N,H)
hold on
loglog(N,(H1Hist),'-*')
loglog(N,(H2Hist),'-*')
loglog(N,(H3Hist),'-*')
loglog(N,(H4Hist),'-*')
loglog(N,(H5Hist),'-*')
loglog(N,(H6Hist),'-*')
% loglog(N,(H7Hist))
title('Hurst exponent with 30db SNR white gaussian noise')
legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
xlabel('Number of Hurst exponents calculated decimated')
ylabel('Hurst exponent')
set(gca,'FontSize', 12);



%% theoretical and estimated hurst exponents with 30dB alpha stable noise
clear all;close all;clc
%hurst parameter
H=0.05:0.05:1;
%fGN generation with variance 1 and 500 samples
sigma=1;
n=4096;
%number of timeseries generated per hurst exponent
N=1:round(100/length(H)):100;    %number of timeseries generated
M=n/2;
force=1;
%Hurst parameter plot calculation
for i=1:length(H)
%     M=20;
    for k=1:length(N)
        auxO=ffGn(n,H(i),1,0);
        aux1 = decimate(auxO,2);
        aux1 = aux1+stblrnd(0.9,1,1,0.1);
     
        H1(k)=hurst_estimate(aux1','RS');
        H2(k)=hurst_estimate(aux1','aggvar');
        H3(k)=hurst_estimate(aux1','boxper');
        H4(k)=hurst_estimate(aux1','diffvar');
        H5(k)=hurst_estimate(aux1','absval');
        H6(k)=hurst_estimate(aux1','peng');
    end
    H1Hist(i)=mean(H1);
    H2Hist(i)=mean(H2);
    H3Hist(i)=mean(H3);
    H4Hist(i)=mean(H4);
    H5Hist(i)=mean(H5);
    H6Hist(i)=mean(H6);
    
end

Horig=H';
H_RS=H1Hist'; H_aggvar=H2Hist'; H_boxper=H3Hist';
H_diffvar=H4Hist'; H_absval=H5Hist';    H_peng=H6Hist';
T = table(Horig, H_RS, H_aggvar, H_boxper, H_diffvar, H_absval,H_peng)

Horig=H';
H_RS=H1Hist'; H_aggvar=H2Hist'; H_boxper=H3Hist';
H_diffvar=H4Hist'; H_absval=H5Hist';    H_peng=H6Hist';
T = table(Horig, H_RS, H_aggvar, H_boxper, H_diffvar, H_absval,H_peng);

%estimation error calculation
TNew=table2array(T);
TNew(isnan(TNew))=0;
[row col]=size(TNew);
TError=TNew;
for j=2:col
   TError(:,j)=abs(((TNew(:,j)-TNew(:,1)))./TNew(:,1))*100;
end
TError
figure()
for i=2:col
    plot(TError(:,1),TError(:,i))
    legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
    hold on
end

figure()
subplot(2,1,1)
bar(TError(1:8),TError(1:8,:))
legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
title('Estimation error for 0.05<H<0.4')
xlabel('Hurst exponent')
ylabel('error')
set(gca,'FontSize', 12);

subplot(2,1,2)
bar(TError(9:20),TError(9:end,:))
legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
title('Estimation error for 0.45<H<1')
xlabel('Hurst exponent')
ylabel('error')
set(gca,'FontSize', 12);

%%
% 
% The following figure shows the robustness of the Hurst parameter
% estimators in the presence of symetric alpha stable random noise with 30 dB of
% intensity and $alpha=0.9$. It can be observed that the the presense on
% noise makes the absolute value estimator works better regarding to the
% other estimators.
%
%
Horig=H';
H_RS=H1Hist'; H_aggvar=H2Hist'; H_boxper=H3Hist';
H_diffvar=H4Hist'; H_absval=H5Hist';    H_peng=H6Hist';
T = table(Horig, H_RS, H_aggvar, H_boxper, H_diffvar, H_absval,H_peng)

figure()
plot(N,H)
hold on
loglog(N,(H1Hist),'-*')
loglog(N,(H2Hist),'-*')
loglog(N,(H3Hist),'-*')
loglog(N,(H4Hist),'-*')
loglog(N,(H5Hist),'-*')
loglog(N,(H6Hist),'-*')
% loglog(N,(H7Hist))
legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
title('Hurst exponent with 30db SNR symetric alpha stable noise')
xlabel('Number of Hurst exponents calculated decimated')
ylabel('Hurst exponent')
set(gca,'FontSize', 12);




%% part 3: Reduced sample length H parameter estimation
%%
% 
% In this section jus 2048 points are employed to generate the fGn which is
% different from part 2, where the full length (4096 samples) was generated
% and then subsampled.
%
%  
%% theoretical and estimated hurst exponents on nominal conditions
clc; clear; close all;
H=0.05:0.05:1; %hurst parameter
sigma=1;       %variance
nSelect=1000:1000:4000; %noise samples
figure()
for kapa=1:length(nSelect)
    n=nSelect(kapa);
    %number of timeseries generated per hurst exponent
    N=1:round(100/length(H)):100;    %number of timeseries generated
    %Hurst parameter plot calculation
    for i=1:length(H)
    %     M=20;
        for k=1:length(N)
            aux1=ffGn(n,H(i),1,0);
%             aux1 = awgn(aux1,30);
            aux1 = aux1+stblrnd(0.9,1,1,0.1);
            H1(k)=hurst_estimate(aux1','RS');
            H2(k)=hurst_estimate(aux1','aggvar');
            H3(k)=hurst_estimate(aux1','boxper');
            H4(k)=hurst_estimate(aux1','diffvar');
            H5(k)=hurst_estimate(aux1','absval');
            H6(k)=hurst_estimate(aux1','peng');
        end
        H1Hist(i)=mean(H1); H2Hist(i)=mean(H2);
        H3Hist(i)=mean(H3); H4Hist(i)=mean(H4);
        H5Hist(i)=mean(H5); H6Hist(i)=mean(H6);
    end
    subplot(1,4,kapa)
    plot(N,H)
    hold on
    loglog(N,(H1Hist),'-*'); loglog(N,(H2Hist),'-*')
    loglog(N,(H3Hist),'-*'); loglog(N,(H4Hist),'-*')
    loglog(N,(H5Hist),'-*'); loglog(N,(H6Hist),'-*')
    legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
    titl=strcat('Hurst exponent for n=',string(n));
    title(titl)
    xlabel('Number of Hurst exponents')
    ylabel('Hurst exponent')
    set(gca,'FontSize', 12);
end

%% theoretical and estimated hurst exponents with 30dB gaussian noise
clear all;close all;clc
%hurst parameter
H=0.05:0.05:1;
%fGN generation with variance 1 and 500 samples
sigma=1;
n=2048;
%number of timeseries generated per hurst exponent
N=1:round(100/length(H)):100;    %number of timeseries generated
M=n/2;
force=1;
%Hurst parameter plot calculation
for i=1:length(H)
%     M=20;
    for k=1:length(N)
        aux1=ffGn(n,H(i),1,0);
%         aux1 = awgn(aux1,30);
        aux1 = aux1+stblrnd(0.9,1,1,0.1);
        H1(k)=hurst_estimate(aux1','RS');
        H2(k)=hurst_estimate(aux1','aggvar');
        H3(k)=hurst_estimate(aux1','boxper');
        H4(k)=hurst_estimate(aux1','diffvar');
        H5(k)=hurst_estimate(aux1','absval');
        H6(k)=hurst_estimate(aux1','peng');
    end
    H1Hist(i)=mean(H1);
    H2Hist(i)=mean(H2);
    H3Hist(i)=mean(H3);
    H4Hist(i)=mean(H4);
    H5Hist(i)=mean(H5);
    H6Hist(i)=mean(H6);
    
end
%%
%
%In the presence of random noise with 30dB SNR, the 
%residual covariance and periodogram  methods returns the best hurst
%estimatior performance for the hurst parameter.
%

Horig=H';
H_RS=H1Hist'; H_aggvar=H2Hist'; H_boxper=H3Hist';
H_diffvar=H4Hist'; H_absval=H5Hist';    H_peng=H6Hist';
T = table(Horig, H_RS, H_aggvar, H_boxper, H_diffvar, H_absval,H_peng)


figure()
plot(N,H)
hold on
loglog(N,(H1Hist),'-*')
loglog(N,(H2Hist),'-*')
loglog(N,(H3Hist),'-*')
loglog(N,(H4Hist),'-*')
loglog(N,(H5Hist),'-*')
loglog(N,(H6Hist),'-*')
% loglog(N,(H7Hist))
legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
title('Hurst exponent for multiple Estimators with 30db SNR white gaussian noise')

%% theoretical and estimated hurst exponents with 30dB alpha stable noise
clear all;close all;clc
%hurst parameter
H=0.05:0.05:1;
%fGN generation with variance 1 and 500 samples
sigma=1;
n=2048;
%number of timeseries generated per hurst exponent
N=1:round(100/length(H)):100;    %number of timeseries generated
M=n/2;
force=1;
%Hurst parameter plot calculation
for i=1:length(H)
%     M=20;
    for k=1:length(N)
         aux1=ffGn(n,H(i),1,0);
%         aux1 = decimate(auxO,2);
        aux1 = aux1+stblrnd(0.9,1,1,0.1);
     
        H1(k)=hurst_estimate(aux1','RS');
        H2(k)=hurst_estimate(aux1','aggvar');
        H3(k)=hurst_estimate(aux1','boxper');
        H4(k)=hurst_estimate(aux1','diffvar');
        H5(k)=hurst_estimate(aux1','absval');
        H6(k)=hurst_estimate(aux1','peng');
    end
    H1Hist(i)=mean(H1);
    H2Hist(i)=mean(H2);
    H3Hist(i)=mean(H3);
    H4Hist(i)=mean(H4);
    H5Hist(i)=mean(H5);
    H6Hist(i)=mean(H6);
    
end
%%
% 
% The following figure shows the robustness of the Hurst parameter
% estimators in the presence of symetric alpha stable random noise with 30 dB of
% intensity and $alpha=0.9$. It can be observed that the presense on
% noise makes the aperiodogram and covariance methods works better regarding to the
% other estimators.
%
%
Horig=H';
H_RS=H1Hist'; H_aggvar=H2Hist'; H_boxper=H3Hist';
H_diffvar=H4Hist'; H_absval=H5Hist';    H_peng=H6Hist';
T = table(Horig, H_RS, H_aggvar, H_boxper, H_diffvar, H_absval,H_peng)

figure()
plot(N,H)
hold on
loglog(N,(H1Hist),'-*')
loglog(N,(H2Hist),'-*')
loglog(N,(H3Hist),'-*')
loglog(N,(H4Hist),'-*')
loglog(N,(H5Hist),'-*')
loglog(N,(H6Hist),'-*')
% loglog(N,(H7Hist))
legend('H real','RS','aggvar','boxper','diffvar','absval','peng')
title('Hurst exponent for multiple Estimators with 30db SNR symetric alpha stable noise')

%% Conclusions
% The estimation of the Hurst parameter is not a trivial task and requuires
% considering multiple conditions in order to determine the correct LRD of
% a signal. The best solution is evaluating multiple calculation methods,
% so a global idea about sesionality and LRD can be abstracted.
%
% Also, the covariance method shows to be a good estimator under most of
% the evaluated conditions in this work like random noise (gaussian and
% alpha stable), decimation, and less sample data. However, the assesment
% of multiple calculation methods is required in order to obtain the right
% estimation and a better knowledge of the system.
%  
%  
% 
