function [] = plotTimeseriesAtLevel(varname, level, plotAll)

y = load([varname,'.dat']);
load('time.dat')

figure;
if time(end) > 3 && ~plotAll
  i = time > 3;
else
  i = 1:length(time);
end
plot(time(i), y(i,level),'linewidth', 2.0)
xlabel('Time [d]', 'fontsize',15)
title(varname,'fontsize',15)


end