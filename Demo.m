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
SI = zeros(num, 3);
Gini = zeros(num, 3);

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
    
    SI(i, 1) = prod(power(SE1, 1/Len)) / mean(SE1);
    SI(i, 2) = prod(power(SE2, 1/Len)) / mean(SE2);
    
    [ASE1, ~] = sort(SE1);
    value1 = 0;
    for k = 1 : Len
        value1 = value1 + 2 * ASE1(k) / norm(ASE1, 1) * (Len - k + 0.5) / Len;
    end
    Gini(i, 1) = 1 - value1;
    
    [ASE2, ~] = sort(SE2);
    value2 = 0;
    for k = 1 : Len
        value2 = value2 + 2 * ASE2(k) / norm(ASE2, 1) * (Len - k + 0.5) / Len;
    end
    Gini(i, 2) = 1 - value2;
    
    disp(['i = ',num2str(i)])
end

figure;
subplot(231);
plot(Kurtosis(:, 1), '-o'); hold on; 
plot(Kurtosis(:, 2), '-d');
xlabel('Amplitude'); ylabel('Value');
ylim([-100, 2000]); 
title('Kurtosis');
legend('TI', 'RTI');
set(gca,'FontName','Times new roman');
subplot(232);
plot(L2L1(:, 1), '-o'); hold on; 
plot(L2L1(:, 2), '-d');
xlabel('Amplitude'); ylabel('Value');
title('L2/L1 norm');
legend('TI', 'RTI');
set(gca,'FontName','Times new roman');
subplot(233); 
plot(NE(:, 1), '-o'); hold on; 
plot(NE(:, 2), '-d');
xlabel('Amplitude'); ylabel('Value');
title('Negative entropy');
legend('TI', 'RTI');
set(gca,'FontName','Times new roman');
subplot(234);
plot(-SI(:, 1), '-o'); hold on; 
plot(-SI(:, 2), '-d');
xlabel('Amplitude'); ylabel('Value');
title('Smoothness index');
legend('TI', 'RTI');
set(gca,'FontName','Times new roman');
subplot(235);
plot(Gini(:, 1), '-o'); hold on; 
plot(Gini(:, 2), '-d');
xlabel('Amplitude'); ylabel('Value');
title('Gini index');
legend('TI', 'RTI');
set(gca,'FontName','Times new roman');
subplot(236);
plot(Res(:, 1), '-o'); hold on; 
plot(Res(:, 2), '-d');
ylim([-0.05, 0.35]); title('Proposed HI');
legend('TI', 'RTI');
xlabel('Amplitude'); ylabel('Value');
set(gca,'FontName','Times new roman');
set(gcf,'Position',[500,300,750,350]);

