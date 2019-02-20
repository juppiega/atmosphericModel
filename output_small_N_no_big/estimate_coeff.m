function [] = estimate_coeff()

load('time.dat')
load('RH.dat')
load('h.dat')
[LCL, top] = find_LCL_top(RH);

load('N_drops.dat')
load('LWC.dat')
load('drop_area.dat')
ind = time > 3;
N_drops = N_drops(ind,:);
LCL = LCL(ind);
top = top(ind);
LWC = LWC(ind,:);
drop_area = drop_area(ind,:);
time = time(ind);
N = length(time);

N_mean = zeros(N,1);
LWP = zeros(N,1);
tau = zeros(N,1);
dz = h(2) - h(1);

for i = 1:N
    N_mean(i) = mean(N_drops(i,LCL(i):top(i)));
    LWP(i) = dz*sum(LWC(i,LCL(i):top(i)));
    tau(i) = dz*2*sum(drop_area(i,LCL(i):top(i)));
end

tau_pred = N_mean.^(1/3) .* LWP.^(5/6);
tau_pred = tau_pred * mean(tau./tau_pred);

[c, c_int] = regress(log(tau), [ones(N,1), log(N_mean), log(LWP)]);
a = c(2)
b = c(3)
r = b/a
a_sig = (c_int(2,2) - a)/2;
b_sig = (c_int(3,2) - b)/2;
r_sig = sqrt((b_sig/a)^2 + (b*a_sig/a^2)^2)
tau_best = exp(c(1)) * N_mean.^a .* LWP.^b;

figure
plot(time, tau, time, tau_pred, time, tau_best)
legend('tau', 'tau monodisp.', 'tau best fit')
corrcoef(tau, tau_best)

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