clear; % Script to generate a plot similar to Fig2 in the Schultz & Telesca (2017) paper.

% Load data and predefine.
T=1:500;
load('CleanTraces.mat');

% Contaminate traces with Gaussian noise and correlate.
sig=0.15;
x15a=X1+normrnd(0, sig, size(T)); x15b=X2+normrnd(0, sig, size(T));
[Tx, X15, C15]=ReshuffleCorr(T,x15b,x15a, 10000,[0.9545,0.9973]);
fprintf('%f: (%f %f)\n',sig, X15(Tx==lag), X15(Tx==lag)/C15(Tx==lag,1));

sig=0.20;
x20a=X1+normrnd(0, sig, size(T)); x20b=X2+normrnd(0, sig, size(T));
[~, X20, C20]=ReshuffleCorr(T,x20b,x20a, 10000,[0.9545,0.9973]);
fprintf('%f: (%f %f)\n',sig, X20(Tx==lag), X20(Tx==lag)/C20(Tx==lag,1));

sig=0.25;
x25a=X1+normrnd(0, sig, size(T)); x25b=X2+normrnd(0, sig, size(T));
[~, X25, C25]=ReshuffleCorr(T,x25b,x25a, 10000,[0.9545,0.9973]);
fprintf('%f: (%f %f)\n',sig, X25(Tx==lag), X25(Tx==lag)/C25(Tx==lag,1));

sig=0.30;
x30a=X1+normrnd(0, sig, size(T)); x30b=X2+normrnd(0, sig, size(T));
[~, X30, C30]=ReshuffleCorr(T,x30b,x30a, 10000,[0.9545,0.9973]);
fprintf('%f: (%f %f)\n',sig, X30(Tx==lag), X30(Tx==lag)/C30(Tx==lag,1));

% Plot.
figure(1); clf;
PURP=[133,57,227]/256;
ms=2;

% s15 plots.
subplot(4,4,[1 2]);
plot(T, x15a+1, '-o', 'MarkerSize',ms); hold on;
plot(T, x15b, '-o', 'MarkerSize',ms);
xlabel('Samples'); ylabel('Ampitude');

subplot(4,4,[3 4]);
plot(Tx, X15, 'Color', PURP); hold on;
plot(Tx, C15(:,1),'-','Color','k');
plot(Tx, C15(:,2),'-','Color','k');
xlabel('Sample Lag'); ylabel('Correlation Ampitude');

% s20 plots.
subplot(4,4,[5 6]);
plot(T, x20a+1, '-o', 'MarkerSize',ms); hold on;
plot(T, x20b, '-o', 'MarkerSize',ms);
xlabel('Samples'); ylabel('Ampitude');

subplot(4,4,[7 8]);
plot(Tx, X20, 'Color', PURP); hold on;
plot(Tx, C20(:,1),'-','Color','k');
plot(Tx, C20(:,2),'-','Color','k');
xlabel('Sample Lag'); ylabel('Correlation Ampitude');

% s25 plots.
subplot(4,4,[9 10]);
plot(T, x25a+1, '-o', 'MarkerSize',ms); hold on;
plot(T, x25b, '-o', 'MarkerSize',ms);
xlabel('Samples'); ylabel('Ampitude');

subplot(4,4,[11 12]);
plot(Tx, X25, 'Color', PURP); hold on;
plot(Tx, C25(:,1),'-','Color','k');
plot(Tx, C25(:,2),'-','Color','k');
xlabel('Sample Lag'); ylabel('Correlation Ampitude');

% s30 plots.
subplot(4,4,[13 14]);
plot(T, x30a+1, '-o', 'MarkerSize',ms); hold on;
plot(T, x30b, '-o', 'MarkerSize',ms);
xlabel('Samples'); ylabel('Ampitude');

subplot(4,4,[15 16]);
plot(Tx, X30, 'Color', PURP); hold on;
plot(Tx, C30(:,1),'-','Color','k');
plot(Tx, C30(:,2),'-','Color','k');
xlabel('Sample Lag'); ylabel('Correlation Ampitude');







