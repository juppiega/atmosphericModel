function plotAtTime(varname, days, hours, minutes)
  % varname = ['u', 'v', 'theta'], time is difference from 0.
  
  load('h.dat')
  if strcmpi(varname,'u')
    load('ua.dat')
    var = ua;
  elseif strcmpi(varname,'v')
    load('va.dat')
    var = va;
  elseif strcmpi(varname,'theta')
    load('theta.dat')
    var = theta-273.15;
  elseif strcmpi(varname,'alpha_pinene')
    load('alpha_pinene.dat')
    var = alpha_pinene;
  elseif strcmpi(varname,'isoprene')
    load('isoprene.dat')
    var = isoprene;
  elseif strcmpi(varname,'Km')
    load('Km.dat')
    var = Km;
    h = (h(2:end) + h(1:end-1)) / 2;
  elseif strcmpi(varname,'Kh')
    load('Kh.dat')
    var = Kh;
    h = (h(2:end) + h(1:end-1)) / 2;
  elseif strcmpi(varname,'Ri')
    load('Ri.dat')
    var = Ri;
    h = (h(2:end) + h(1:end-1)) / 2;
  else
    error('Unrecognized variable name')
  end
  
  load('time.dat')
  
  
  plotTime = days + hours/24 + minutes/60;
  if plotTime > max(time) || plotTime < min(time)
    error(['Plot time not in the domain of model output times: [', num2str(time(1)),', ', num2str(time(end)),'] days'])
  end
  [~,plotInd] = min(abs(time - plotTime));
  
  varAtTime = var(plotInd,:);
  
  figure;
  plot(varAtTime, h, 'linewidth', 2.0)
  
  title([varname,' at ',num2str(days),' d ', num2str(hours),' h ', num2str(minutes),' m'], 'fontsize', 15)
  xlabel(varname, 'fontsize', 15)
  ylabel('z [m]', 'fontsize', 15)
  ylim([0, 3000])
  set(gca,'fontsize',15)
  
end