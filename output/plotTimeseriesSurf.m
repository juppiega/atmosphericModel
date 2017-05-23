function plotTimeseriesSurf(varname, chemname)

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
  elseif strcmpi(varname,'chemistry')
    var = load([chemname,'.dat']);
    var = var(time >= 3, :);
    time = time(time > 3);
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
  if nargin == 1
    title(varname, 'fontsize', 15)
  else
    title(chemname, 'fontsize', 15)
  end
  set(gca,'fontsize',15)
  ylim([0, 3000])
  
  
  

end