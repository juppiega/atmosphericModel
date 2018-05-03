function plotTimeseriesSurf(varname, chemname)

  load('h.dat')
  load('time.dat')
  size_distrib = false;
  if strcmpi(varname,'u')
    load('ua.dat')
    var = ua;
  elseif strcmpi(varname,'v')
    load('va.dat')
    var = va;
  elseif strcmpi(varname,'q')
    load('q.dat')
    var = q;
  elseif strcmpi(varname,'theta')
    load('theta.dat')
    var = theta;
  elseif strcmpi(varname,'chemistry')
    var = load([chemname,'.dat']);
    if strcmpi(chemname, 'size_distribution')
      var = log10(var/1E6);
      var(var < 0) = 0;
    end  
    var = var(time >= 3, :);
    time = time(time >= 3);
  elseif strcmpi(varname,'Km')
    load('Km.dat')
    var = Km;
    h = (h(2:end) + h(1:end-1)) / 2;
    time = time(2:end);
  elseif strcmpi(varname,'Kh')
    load('Kh.dat')
    var = (Kh);
    h = (h(2:end) + h(1:end-1)) / 2;
    time = time(2:end);
  elseif strcmpi(varname,'Ri')
    load('Ri.dat')
    var = log10(abs(Ri));
    h = (h(2:end) + h(1:end-1)) / 2;
    time = time(2:end);
    varname = 'log10(abs(Ri))'
  elseif strcmpi(varname,'t_diff')
    theta = load('theta.dat');
    theta_u = load('theta_u.dat');
    var = theta_u - theta(2:end,1:end-1);
    var(var < -200) = 0;
    h = (h(2:end) + h(1:end-1)) / 2;
    time = time(2:end);
  else
    var = load([varname,'.dat']);
    h = (h(2:end) + h(1:end-1)) / 2;
    time = time(2:end);
  end
  
  if length(time) < size(var,1)
    time = [0; time];
  end
  if length(h) ~= size(var,2)
    diameter(1)=2D-9;
    nr_bins = size(var,2);
    for i=2:nr_bins
      diameter(i)=diameter(i-1)*(2.5D-6/diameter(1))**(1D0/(nr_bins-1));
    end
    h = diameter * 1E9;
    size_distrib = true;
  end
 
  
  [X,Y] = meshgrid(time, h);
  
  figure;
  surf(X, Y, var','edgecolor','none')
  view(2)
  axis tight
  colorbar
  if strcmpi(varname,'w_kin') || strcmpi(varname,'flux_t')
      caxis([-max(var(:)), max(var(:))])
  end
  if size_distrib
    set(gca, 'YScale', 'log')
  end
  xlabel('Time [d]', 'fontsize', 15)
  if size_distrib
    ylabel('Diameter [nm]', 'fontsize', 15)
  else
    ylabel('z [m]', 'fontsize', 15)
  end
  if nargin == 1
    title(varname, 'fontsize', 15)
  else
    title(chemname, 'fontsize', 15)
  end
  set(gca,'fontsize',15)
  if ~size_distrib
    ylim([0, 3000])
  end
  
  
  

end