clear; clc; close all;

%% Parâmetros 
Kp   = 0.8;      % ganho do integrador K_P* [1/s]
m0   = 42;       % valor inicial da entrada
c0   = 60;       % valor inicial da saída
dm   = 5;        % amplitude do pulso (+5 e depois -5)
Ts   = 0.01;     % passo de tempo [s]
t    = (0:Ts:20)';

%% Entrada
m = m0 * ones(size(t));
m(t >= 5  & t < 10) = m0 + dm;
m(t >= 10 & t < 15) = m0 - dm;
m1 = m0 + dm;    % amplitude do 1º pulso
m2 = m0 - dm;    % amplitude do 2º pulso


%% CASO 1 -INTEGRADOR PURO

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

%%  CASO 2 - INTEGRADOR REAL 
tau   = 0.5;     
theta = 0.8;     

Gp  = tf(Kp, [tau 1 0], 'InputDelay', theta);   
Gd  = tf(Kp, [tau 1],   'InputDelay', theta);  

c2  = c0 + lsim(Gp, m - m0, t);        % saída
dc2 = lsim(Gd, m - m0, t);             % derivada da saída

plotFig(t, m, c2, dc2, [-4 4], 'east');
sgtitle('Integrador real');


s1r = interp1(t, dc2, 10 + theta - Ts);   % fim do 1º pulso (já com atraso)
s2r = interp1(t, dc2, 15 + theta - Ts);   % fim do 2º pulso (já com atraso)
Kp_est_r = (s2r - s1r) / (m2 - m1);

fprintf('--- Integrador real (tau = %.2f s, theta = %.2f s) ---\n', tau, theta);
fprintf('s1 = %.3f | s2 = %.3f | m1 = %.1f | m2 = %.1f\n', s1r, s2r, m1, m2);
fprintf('Kp* estimado = %.3f (valor real: %.3f)\n', Kp_est_r, Kp);


%% IDENTIFICAÇAO PELO MÉTODO FIT3 

t_on  = 5;                                
t_off = 10;                               

dcm = gradient(c2, Ts);                   % derivada da saída medida c'(t)

% Dm e Dcs
Dm  = m1 - m0;                                    % degrau na entrada
dc0 = mean(dcm(t < t_on));                        % patamar inicial de c'
dcf = dcm(find(t >= t_off, 1) - 1);               % patamar de c' antes da descida
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

%% Modelo identificado x dados medidos

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


%%  SINTONIA DO CONTROLADOR PROPOSTO E MALHA FECHADA lambda

lambda = t0_id;                           % lambda = t0
Kc_l   = 1 / (K_id * (lambda + t0_id));   % Kc'
tauD_l = tau_id;                          % tauD'
alpha  = 5;                               % coeficiente do filtro derivativo

Kc_PI = tau_id / (K_id * (lambda + t0_id));  %DÚVIDA: a equação (13) do artigo tá diferente                  % Kc da forma PI 
tauI  = tau_id;                                                % tau_i da forma PI
C_PD  = tf(Kc_l * [(alpha+1)*tauD_l, 1], [alpha*tauD_l, 1]);   % C(s) 
%DÚVIDA: não entendi de onde vem o termo "(alpha+1)*tauD_l, 1"

fprintf('\n================ CONTROLADOR PROPOSTO ================\n');
fprintf('Modelo identificado:\n');
fprintf('  G*_P(s) = %.4f e^(-%.4f s) / (%.4f s + 1)\n', K_id, t0_id, tau_id);
fprintf('  lambda = t0 = %.4f s | alpha = %g\n\n', lambda, alpha);

fprintf('Forma PI x s (Eq. 13):   Gc(s) = Kc (1 + 1/(tau_i s)) s\n');
fprintf('  Kc    = tau / (K (lambda + t0)) = %.4f\n', Kc_PI);
fprintf('  tau_i = tau                     = %.4f s\n\n', tauI);

fprintf('Forma PD (Eqs. 14 e 15): C(s) = Kc'' (1 + tauD'' s / (alpha tauD'' s + 1))\n');
fprintf('  Kc''   = Kc / tau_i = 1 / (K (lambda + t0)) = %.4f\n', Kc_l);
fprintf('  tauD'' = tau_i = tau                        = %.4f s\n\n', tauD_l);

fprintf('Função de transferência do controlador:\n');
fprintf('  C(s) = (%.4f s + %.4f) / (%.4f s + 1)\n\n', ...
        Kc_l*(alpha+1)*tauD_l, Kc_l, alpha*tauD_l);
disp(C_PD);
fprintf('======================================================\n\n');


c_sp = 70;                                % valor final do setpoint
t_sp = 2;                                 % instante do degrau no setpoint [s]
r    = c0 + (c_sp - c0) * (t >= t_sp);    % setpoint

N   = numel(t);
nd  = round(theta / Ts);                  % tempo morto em amostras
a   = exp(-Ts / tau);                     % discretização da inércia real

y   = zeros(N,1);  y(1) = c0;             % variável controlada
ud  = zeros(N,1);                         % desvio da variável manipulada
x   = 0;                                  % estado: derivada da saída do processo
z   = 0;                                  % estado do filtro do termo derivativo

for k = 1:N
    e     = r(k) - y(k);
    D     = (e - z) / alpha;              % = tauD s/(alpha tauD s + 1) aplicado a e
    ud(k) = Kc_l * (e + D);               % ação de controle 
    z     = z + Ts * (e - z) / (alpha * tauD_l);

    if k > nd, uin = ud(k - nd); else, uin = 0; end   % tempo morto do processo
    x = a * x + (1 - a) * Kp * uin;       % inércia do processo (dy/dt)
    if k < N, y(k+1) = y(k) + Ts * x; end % integrador do processo
end
u = m0 + ud;                              % variável manipulada

%% Índice IAE
e_tot   = r - y;
IAE_sum = sum(abs(e_tot));                % soma das amostras
IAE_int = trapz(t, abs(e_tot));           % integral em [s]
fprintf('--- Método proposto ---\n');
fprintf('IAE = %.4e (soma das amostras) | %.4f (integral)\n', IAE_sum, IAE_int);


figure('Color','w','Position',[100 100 900 380]);
tiledlayout(1, 2, 'TileSpacing','compact', 'Padding','compact');
lw = 1.5;

nexttile;
plot(t, r, 'LineWidth', lw); hold on;
plot(t, y, 'LineWidth', lw);
title('(a)', 'FontWeight','bold');
xlabel('Time [s]'); ylabel('Controlled Variable');
legend('Setpoint', 'Proposed', 'Location','southeast');
xlim([0 20]); ylim([58 74]); xticks(0:5:20); yticks(60:2:74); box on;

nexttile;
plot(t, u, 'LineWidth', lw);
title('(b)', 'FontWeight','bold');
xlabel('Time [s]'); ylabel('Manipulated Variable');
legend('Proposed', 'Location','northeast');
xlim([0 20]); ylim([36 54]); xticks(0:5:20); yticks(38:2:52); box on;



function plotFig(t, m, c, dc, ylimC, legLocC)
    figure('Color','w','Position',[100 100 900 550]);
    tiledlayout(2, 2, 'TileSpacing','compact', 'Padding','compact');
    lw = 1.5;

    % (a) Saída
    nexttile(1);
    plot(t, c, 'LineWidth', lw);
    title('(a)', 'FontWeight','bold');
    xlabel('Time [s]'); ylabel('Amplitude');
    legend('Output, c(t)', 'Location','northeast');
    xlim([0 20]); ylim([60 80]); xticks(0:5:20); box on;

    % (b) Entrada
    nexttile(3);
    plot(t, m, 'LineWidth', lw);
    title('(b)', 'FontWeight','bold');
    xlabel('Time [s]'); ylabel('Amplitude');
    legend('Input, m(t)', 'Location','southeast');
    xlim([0 20]); ylim([36 47.5]); xticks(0:5:20); yticks(38:2:46); box on;

    % (c) Derivada da saída (ocupa a coluna da direita inteira)
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