clear; clc; close all;

Kp   = 0.8;
m0   = 42;
c0   = 60;
dm   = 5;
Ts   = 0.01;
t    = (0:Ts:20)';

m = m0 * ones(size(t));
m(t >= 5  & t < 10) = m0 + dm;
m(t >= 10 & t < 15) = m0 - dm;
m1 = m0 + dm;
m2 = m0 - dm;

c1  = c0 + Kp * cumtrapz(t, m - m0);
dc1 = Kp * (m - m0);

plotFig(t, m, c1, dc1, [-5 5], 'southeast');
sgtitle('Integrador puro');

i1 = t >= 5  & t < 10;   p1 = polyfit(t(i1), c1(i1), 1);   s1 = p1(1);
i2 = t >= 10 & t < 15;   p2 = polyfit(t(i2), c1(i2), 1);   s2 = p2(1);
Kp_est = (s2 - s1) / (m2 - m1);

fprintf('--- Integrador puro ---\n');
fprintf('s1 = %.3f | s2 = %.3f | m1 = %.1f | m2 = %.1f\n', s1, s2, m1, m2);
fprintf('Kp* estimado = %.3f (valor real: %.3f)\n\n', Kp_est, Kp);

tau   = 0.5;
theta = 0.8;

Gp  = tf(Kp, [tau 1 0], 'InputDelay', theta);
Gd  = tf(Kp, [tau 1],   'InputDelay', theta);

c2  = c0 + lsim(Gp, m - m0, t);
dc2 = lsim(Gd, m - m0, t);

plotFig(t, m, c2, dc2, [-4 4], 'east');
sgtitle('Integrador real');

s1r = interp1(t, dc2, 10 + theta - Ts);
s2r = interp1(t, dc2, 15 + theta - Ts);
Kp_est_r = (s2r - s1r) / (m2 - m1);

fprintf('--- Integrador real (tau = %.2f s, theta = %.2f s) ---\n', tau, theta);
fprintf('s1 = %.3f | s2 = %.3f | m1 = %.1f | m2 = %.1f\n', s1r, s2r, m1, m2);
fprintf('Kp* estimado = %.3f (valor real: %.3f)\n', Kp_est_r, Kp);

t_on  = 5;
t_off = 10;

dcm = gradient(c2, Ts);

Dm  = m1 - m0;
dc0 = mean(dcm(t < t_on));
dcf = dcm(find(t >= t_off, 1) - 1);
Dcs = dcf - dc0;

ok  = t < t_off;
t1  = tcross(t(ok), dcm(ok), dc0 + 0.283*Dcs) - t_on;
t2  = tcross(t(ok), dcm(ok), dc0 + 0.632*Dcs) - t_on;

K_id   = Dcs / Dm;
tau_id = 1.5 * (t2 - t1);
t0_id  = t2 - tau_id;

fprintf('\n--- Identificação Fit3 (a partir do pulso) ---\n');
fprintf('Dm = %.2f | Dcs = %.4f | t1 = %.4f s | t2 = %.4f s\n', Dm, Dcs, t1, t2);
fprintf('K = %.4f | tau = %.4f s | t0 = %.4f s\n', K_id, tau_id, t0_id);
fprintf('G*_P(s) = %.4f e^(-%.4f s) / (%.4f s + 1)\n', K_id, t0_id, tau_id);
fprintf('G_P(s)  = %.4f e^(-%.4f s) / (s(%.4f s + 1))\n', K_id, t0_id, tau_id);

Gstar_id = tf(K_id, [tau_id 1], 'InputDelay', t0_id);
dc_fit   = dc0 + lsim(Gstar_id, m - m0, t);

figure('Color','w','Position',[100 100 700 450]);
plot(t, dcm, 'LineWidth', 1.8); hold on;
plot(t, dc_fit, '--r', 'LineWidth', 1.5);
yline(dc0 + 0.283*Dcs, ':k'); yline(dc0 + 0.632*Dcs, ':k');
xline(t_on + t1, ':k');       xline(t_on + t2, ':k');
plot([t_on+t1 t_on+t2], [dc0+0.283*Dcs dc0+0.632*Dcs], 'ko', 'MarkerFaceColor','k');
text(t_on+t1, dc0+0.283*Dcs, '  0,283\Delta c_s', 'VerticalAlignment','top');
text(t_on+t2, dc0+0.632*Dcs, '  0,632\Delta c_s', 'VerticalAlignment','top');
xlabel('Time [s]'); ylabel('Amplitude'); box on;
xlim([0 20]);
legend('c''(t) medida', 'FOPDT identificado (Fit3)', 'Location','southeast');
title(sprintf('Fit3: K = %.4f, \\tau = %.4f s, t_0 = %.4f s', K_id, tau_id, t0_id));

rng(0);
sigma_n = 0.12;
tau_d   = 1.0;

c_ruido  = c2 + sigma_n * randn(size(t));

Gder     = tf([1 0], [tau_d 1]);
dc_ruido = lsim(Gder, c_ruido - c0, t);

plotFig(t, m, c_ruido, dc_ruido, [-5 5], 'southwest');
sgtitle('Integrador real com ruído de medição');

function plotFig(t, m, c, dc, ylimC, legLocC)
    figure('Color','w','Position',[100 100 900 550]);
    tiledlayout(2, 2, 'TileSpacing','compact', 'Padding','compact');
    lw = 1.5;

    nexttile(1);
    plot(t, c, 'LineWidth', lw);
    title('(a)', 'FontWeight','bold');
    xlabel('Time [s]'); ylabel('Amplitude');
    legend('Output, c(t)', 'Location','northeast');
    xlim([0 20]); ylim([60 80]); xticks(0:5:20); box on;

    nexttile(3);
    plot(t, m, 'LineWidth', lw);
    title('(b)', 'FontWeight','bold');
    xlabel('Time [s]'); ylabel('Amplitude');
    legend('Input, m(t)', 'Location','southeast');
    xlim([0 20]); ylim([36 47.5]); xticks(0:5:20); yticks(38:2:46); box on;

    nexttile(2, [2 1]);
    plot(t, dc, 'LineWidth', lw);
    title('(c)', 'FontWeight','bold');
    xlabel('Time [s]'); ylabel('Amplitude');
    legend('Output Derivative, c''(t)', 'Location', legLocC);
    xlim([0 20]); ylim(ylimC); xticks(0:5:20); yticks(ylimC(1):1:ylimC(2)); box on;
end

function tc = tcross(t, y, nivel)
    i = find(y >= nivel, 1, 'first');
    tc = t(i-1) + (nivel - y(i-1)) * (t(i) - t(i-1)) / (y(i) - y(i-1));
end