clear;
close all;
clc;

load NASA_IMS_1_7.mat;
load NASA_HI_1_7.mat;
Data = NASA_IMS_1_7;
[Len, col] = size(Data);

Kurtosis = zeros(1, col);
RMS = zeros(1, col);
L2L1 = zeros(1, col);
SI = zeros(1, col);
Gini = zeros(1, col);

for i = 1 : col
    
    data = Data(:, i);
    
    Kurtosis(i) = kurtosis(data);
    
    RMS(i) = rms(data);
    
    L2L1(i) = norm(data, 2) / norm(data, 1);
    
    SE = (abs(hilbert(data))).^2;
    
    SI(i) = prod(power(SE, 1/Len)) / mean(SE);
    
    [ASE, ~] = sort(SE);
    value = 0;
    for k = 1 : Len
        value = value + 2 * ASE(k) / norm(ASE, 1) * (Len - k + 0.5) / Len;
    end
    Gini(i) = 1 - value;
    
    disp(['i = ',num2str(i)])
end
index = 1 : col;
figure;
subplot(231);
plot(index, Kurtosis, 'b');
xlim([1, col]); ylim([0, 60]);
xlabel('Time (×10min)');
ylabel('Value');
title('Kurtosis');
set(gca,'FontName','Times new roman');
subplot(232);
plot(index, RMS, 'b');
xlim([1, col]); ylim([0.1, 0.22]);
xlabel('Time (×10min)');
ylabel('Value');
title('RMS');
set(gca,'FontName','Times new roman');
subplot(233);
plot(index, L2L1, 'b');
xlim([1, col]); ylim([0.0075, 0.0097]);
xlabel('Time (×10min)');
ylabel('Value');
title('L2/L1 norm');
set(gca,'FontName','Times new roman');
subplot(234);
plot(index, -SI, 'b');
xlim([1, col]); ylim([-0.7, -0.35]);
xlabel('Time (×10min)');
ylabel('Value');
title('Smoothness index');
set(gca,'FontName','Times new roman');
subplot(235);
plot(index, Gini, 'b');
xlim([1, col]); ylim([0.35, 0.7]);
xlabel('Time (×10min)');
ylabel('Value');
title('Gini index');
set(gca,'FontName','Times new roman');
subplot(236);
plot(index, res, 'b');
xlim([1, col]); ylim([-0.02, 0.65]);
xlabel('Time (×10min)');
ylabel('Value');
title('Proposed HI');
set(gca,'FontName','Times new roman');
set(gcf,'Position',[500,300,800,350]);




















