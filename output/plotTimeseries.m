function plotTimeseries(varname, beginDay, endDay)
% varname = [chem_param]

load('time.dat')
if strcmpi(varname, 'chem_param')
    alpha_pinene = load('alp_em.dat');
    isoprene = load('is_em.dat');
end

ind = beginDay <= time & time <= endDay;
chemicalComponents = [alpha_pinene, isoprene];

figure;
plotTime = repmat(time(ind), 1, size(chemicalComponents,2));
plot(plotTime, chemicalComponents(ind,:),'linewidth',2.0)
legend('alpha pinene', 'isoprene')
set(gca,'fontsize', 15)
xlabel('Time [days]', 'fontsize', 15)

end