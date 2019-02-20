function [] = compute_indirect()

load('output_small_N/time.dat')
load('output_small_N/RH.dat')
load('output_small_N/h.dat')
[LCL_small, top_small] = find_LCL_top(RH);

load('output_large_N/RH.dat')
[LCL_large, top_large] = find_LCL_top(RH);

ind = time >= 2;
LCL_small = LCL_small(ind);
LCL_large = LCL_large(ind);
top_small = top_small(ind);
top_large = top_large(ind);

load('output_small_N/N_drops.dat')
load('output_small_N/r_eff.dat')
load('output_small_N/LWC.dat')
load('output_small_N/drop_area.dat')
load('output_small_N/rain_rate.dat')
load('output_small_N/q.dat')
load('output_small_N/condensation.dat')
N_small = N_drops(ind,:);
r_eff_small = r_eff(ind,:);
LWC_small = LWC(ind,:);
drop_area_small = drop_area(ind,:);
rain_rate_small = rain_rate(ind,:);
q_small = q(ind,:);
cond_small = condensation(ind,:);

load('output_large_N/N_drops.dat')
load('output_large_N/r_eff.dat')
load('output_large_N/LWC.dat')
load('output_large_N/drop_area.dat')
load('output_large_N/rain_rate.dat')
load('output_large_N/q.dat')
load('output_large_N/condensation.dat')
N_large = N_drops(ind,:);
r_eff_large = r_eff(ind,:);
LWC_large = LWC(ind,:);
drop_area_large = drop_area(ind,:);
rain_rate_large = rain_rate(ind,:);
q_large = q(ind,:);
cond_large = condensation(ind,:);

time = time(ind);

[N_mean_small, LWP_small, tau_small, rain_small, q_int_small, r_mean_small, cond_int_small] = ...
    compute_N_LWP(N_small, LWC_small, drop_area_small, rain_rate_small, q_small, r_eff_small, cond_small, LCL_small, top_small, h);
[N_mean_large, LWP_large, tau_large, rain_large, q_int_large, r_mean_large, cond_int_large] = ...
    compute_N_LWP(N_large, LWC_large, drop_area_large, rain_rate_large, q_large, r_eff_large, cond_large, LCL_large, top_large, h);

dLWP_lar = [0;diff(LWP_large)] / (time(2)-time(1)); % mm/day
dLWP_sma = [0;diff(LWP_small)] / (time(2)-time(1));

conv_start = 0.3; % Fractional days after midnight

RH = RH(ind,:);
LWC_change = onlyCloud(LWC_large - LWC_small, RH);
figure;
loglog(abs(LWC_change(:)), abs(N_small(:)),'.');
r=corr(abs(LWC_change(:)), abs(N_small(:)), 'type','spearman');
legend(['r_{rank} = ', num2str(r)])
xlabel('LWC change')
ylabel('N')

i = time - time(1) > 0.3;
dt = time(2)-time(1);
time = (time-time(1))*24;
figure;
set(gcf,'color','w');
subplot(2,2,1)
plot(time, r_mean_large*1e6, time, r_mean_small*1e6, 'linewidth', 2.0);
title('r_{eff} [\mu m]')
legend('Upwind', 'PLM','location','northwest')
xlim([0,24])
subplot(2,2,2)
plot(time, LWP_large*1e3, time, LWP_small*1e3, 'linewidth', 2.0);
title('LWP [g/m^2]')
xlim([0,24])
xlabel('time [h]')
%legend('Low CCN', 'High CCN')
subplot(2,2,3)
%plot(time, q_int_small, time, q_int_large);
plot(time, rain_large*86400, time, rain_small*86400, 'linewidth', 2.0);
hold all
plot(get(gca,'xlim'), [0,0], 'k--')
title('Rain @ LCL [mm/day]')
xlim([0,24])
xlabel('time [h]')

rain_diff = dt*sum((rain_large(i) - rain_small(i))*86400)
rain_diff_rel = rain_diff / (dt*sum(rain_small(i)*86400))
% cond_diff = sum(cond_int_large(i) - cond_int_small(i))*86400*dt
% cond_diff_rel = sum(cond_int_large(i) - cond_int_small(i)) / sum(cond_int_small(i))
%legend('Low CCN', 'High CCN')
subplot(2,2,4)
plot(time, tau_large, time, tau_small, 'linewidth', 2.0);
tau_ratio = mean(mean(tau_small./tau_large))
title('\tau')
xlim([0,24])
%legend('Low CCN', 'High CCN')
xlabel('time [h]')

N_day = round(1 / (time(2)-time(1))) + 1;
figure;
r = 1.92;
r_sig = 0.09;
R21 = r*(log(LWP_large) - log(LWP_small)) ./ (log(N_mean_large) - log(N_mean_small));
mean_R21 = mean(R21)
%R21 = mean(reshape(R21, N_day, length(time)/N_day),2);
plot(time, R21, 'linewidth', 2.0)
ylabel('R_{2/1}', 'fontsize', 15)
title('R_{2/1}', 'fontsize', 15)
xlabel('Time [h]', 'fontsize', 15)
set(gca,'fontsize',15)
xlim([0,24])

R21_LWP = corr(R21, log(LWP_large) - log(LWP_small));
%corr(R21, 1./(log(N_mean_large) - log(N_mean_small)))
figure;
plot(R21, log(LWP_large) - log(LWP_small),'.');
legend(['r_{linear} = ', num2str(R21_LWP)],'location','southeast')
xlabel('R_{2/1}','fontsize',15)
ylabel('\Delta ln(LWP)','fontsize',15)
set(gca,'fontsize',15)

figure;
plot(R21,rain_small*86400,'.')
r = corr(R21, rain_small,'type','spearman');
legend(['r_{linear} = ', num2str(r)])

figure;
m=6*4;
r = xcorr(R21, rain_small, m, 'coeff');
plot((-m:m)/4, r)


end

function var = onlyCloud(var, RH)

[LCL, top] = find_LCL_top(RH);

for i = 1:length(LCL)
    v = var(i,LCL(i):top(i));
    var(i,:) = 0;
    var(i,LCL(i):top(i)) = v;
end

end


function [N_mean, LWP, tau, rain, q_int, r_mean, cond_int] = compute_N_LWP(N_drops, LWC, drop_area, rain_rate, q, r_eff, cond, LCL, top, h)

N = length(LCL);
N_mean = zeros(N,1);
LWP = zeros(N,1);
tau = zeros(N,1);
rain = zeros(N,1);
q_int = zeros(N,1);
r_mean = zeros(N,1);
cond_int = zeros(N,1);
dz = h(2) - h(1);

for i = 1:N
    N_mean(i) = mean(N_drops(i,LCL(i):top(i)));
    r_mean(i) = mean(r_eff(i,LCL(i):top(i)));
    LWP(i) = dz*sum(LWC(i,LCL(i):top(i)));
    tau(i) = dz*2*sum(drop_area(i,LCL(i):top(i)));
    rain(i) = rain_rate(i,LCL(i));
    q_int(i) = dz*sum(q(i,LCL(i):top(i)));
    cond_int(i) = dz*sum(cond(i,LCL(i):top(i)));
end

end

function [LCL, top] = find_LCL_top(RH)

N = size(RH,1);
m = size(RH,2);
LCL = zeros(N,1);
top = zeros(N,1);
for i = 1:N
    LCL_this = find(RH(i,2:end) > 1, 1) + 1;
    top_this = find(RH(i,:) > 1, 1, 'last');
    if ~isempty(LCL_this)
        LCL(i) = LCL_this;
    else
        LCL(i) = m;
    end
    if ~isempty(top_this)
        top(i) = top_this;
    else
        top(i) = m;
    end
end

end