clc;
clear;
close all;

% Rotating harmonics with background noise
load 'SN.mat';
% Transient interference denoted by a single impluse
load 'st.mat';
% Recipient transient impluses induced by the local defect
load 'RT.mat';
load 'Res.mat';
time = 1;
% Sampling frequency
fs = 10000;
t= (0:time*fs-1)/fs;
Len = time * fs;

figure;
subplot(311);
plot(t, SN, 'b');
xlim([0, 0.1]); ylim([-2, 2]);
xlabel('Time (s)'); ylabel('Amplitude (g)');
set(gca,'FontName','Times new roman');
subplot(312);
plot(t, st, 'b');
xlim([0, 0.1]); ylim([-2, 2]);
xlabel('Time (s)'); ylabel('Amplitude (g)');
set(gca,'FontName','Times new roman');
subplot(313);
plot(t, RT, 'b');
xlim([0, 0.1]); ylim([-2, 2]);
xlabel('Time (s)'); ylabel('Amplitude (g)');
set(gca,'FontName','Times new roman');
set(gcf,'Position',[500,300,500,400]);

num = 30;
Kurtosis = zeros(num, 3);
NE = zeros(num, 3);
L2L1 = zeros(num, 3);
ApEn = zeros(num, 3);
NM = zeros(num, 3);


for i = 1 : num
    st1 = i * st + SN;
    RT1 = i * RT + SN;
    
    Kurtosis(i, 1) = kurtosis(st1);
    Kurtosis(i, 2) = kurtosis(RT1);
    
    L2L1(i, 1) = norm(st1, 2) / norm(st1, 1);
    L2L1(i, 2) = norm(RT1, 2) / norm(RT1, 1);
    
    SE1 = (abs(hilbert(st1))).^2;
    SE2 = (abs(hilbert(RT1))).^2;
    
    NE(i, 1) = mean(SE1/mean(SE1) .* log(SE1/mean(SE1)));
    NE(i, 2) = mean(SE2/mean(SE2) .* log(SE2/mean(SE2)));
    
%     ApEn(i, 1) = approximateEntropy(st1);
%     ApEn(i, 2) = approximateEntropy(RT1);

    [~, Mob1, ~] = HjorthParameters(st1');
    [~, Mob2, ~] = HjorthParameters(RT1');
    
    NM(i, 1) = -Mob1;
    NM(i, 2) = -Mob2;
    
    disp(['i = ',num2str(i)])
end

figure;
subplot(231);
plot(Kurtosis(:, 1), '-o'); hold on; 
plot(Kurtosis(:, 2), '-d');
xlabel('Amplitude (g)'); ylabel('Kurtosis');
ylim([-100, 2000]); 
leg = legend('TI', 'RTI');
leg.ItemTokenSize = [10,20];
set(gca,'FontName','Times new roman');
subplot(232);
plot(L2L1(:, 1), '-o'); hold on; 
plot(L2L1(:, 2), '-d');
xlabel('Amplitude (g)'); ylabel('L2/L1 norm');
leg = legend('TI', 'RTI');
leg.ItemTokenSize = [10,20];
set(gca,'FontName','Times new roman');
subplot(233); 
plot(NE(:, 1), '-o'); hold on; 
plot(NE(:, 2), '-d');
xlabel('Amplitude (g)'); ylabel('NE');
leg = legend('TI', 'RTI');
leg.ItemTokenSize = [10,20];
set(gca,'FontName','Times new roman');
subplot(234);
load 'ApEn.mat';
plot(-ApEn(:, 1), '-o'); hold on; 
plot(-ApEn(:, 2), '-d');
xlabel('Amplitude (g)'); ylabel('ApEn');
leg = legend('TI', 'RTI');
leg.ItemTokenSize = [10,20];
set(gca,'FontName','Times new roman');
subplot(235);
plot(NM(:, 1), '-o'); hold on; 
plot(NM(:, 2), '-d');
ylim([-1.4, -1.1]);
xlabel('Amplitude (g)'); ylabel('NM');
leg = legend('TI', 'RTI');
leg.ItemTokenSize = [10,20];
set(gca,'FontName','Times new roman');
subplot(236);
plot(Res(:, 1), '-o'); hold on; 
plot(Res(:, 2), '-d');
ylim([-0.05, 0.35]); 
leg = legend('TI', 'RTI');
leg.ItemTokenSize = [10,20];
xlabel('Amplitude (g)'); ylabel('Proposed HI');
set(gca,'FontName','Times new roman');
set(gcf,'Position',[500,300,650,300]);

