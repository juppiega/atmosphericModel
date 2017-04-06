function [] = plotTimeseriesAtLevel(varname, level)

y = load([varname,'.dat']);
load('time.dat')

figure;
plot(time, y(:,level),'linewidth', 2.0)
xlabel('Time [d]', 'fontsize',15)
title(varname,'fontsize',15)


end