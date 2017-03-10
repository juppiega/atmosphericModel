function plotTimeseriesSurf(varname)

  load('h.dat')
  load('time.dat')
  if strcmpi(varname,'u')
    load('ua.dat')
    var = ua;
  elseif strcmpi(varname,'v')
    load('va.dat')
    var = va;
  elseif strcmpi(varname,'theta')
    load('theta.dat')
    var = theta;
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
    time = time(2:end);
  elseif strcmpi(varname,'Kh')
    load('Kh.dat')
    var = Kh;
    h = (h(2:end) + h(1:end-1)) / 2;
    time = time(2:end);
  elseif strcmpi(varname,'Ri')
    load('Ri.dat')
    var = Ri;
    h = (h(2:end) + h(1:end-1)) / 2;
    time = time(2:end);
  else
    error('Unrecognized variable name')
  end
  
  [X,Y] = meshgrid(time, h);
  
  figure;
  surf(X, Y, var','edgecolor','none')
  view(2)
  axis tight
  colorbar
  xlabel('Time [d]', 'fontsize', 15)
  ylabel('z [m]', 'fontsize', 15)
  title(varname, 'fontsize', 15)
  set(gca,'fontsize',15)
  ylim([0, 3000])
  
  
  

end