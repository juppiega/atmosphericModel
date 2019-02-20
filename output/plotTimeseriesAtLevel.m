function [] = plotTimeseriesAtLevel(varname, level)

if strcmpi(varname,'Bowen')
  flux_t = load('flux_t.dat');
  flux_q = load('flux_q.dat');
  y = 1004*flux_t(:,level) ./ (2.5E6*flux_q(:,level));
else
  y = load([varname,'.dat']);
end

load('time.dat')
if size(y,1) ~= length(time)
    time(1) = [];
end

figure;
i = 1:length(time);
plot(time(i), y(i,level),'linewidth', 2.0)
xlabel('Time [d]', 'fontsize',15)
title(varname,'fontsize',15)
set(gca,'fontsize',15)
if strcmpi(varname,'Bowen')
  ylim([-5,10]);
end

end