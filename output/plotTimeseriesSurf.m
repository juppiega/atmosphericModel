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
  elseif strcmpi(varname,'PN')
    load('PN.dat')
    var = PN;
  elseif strcmpi(varname,'PM')
    load('PM.dat')
    var = PM;
  elseif strcmpi(varname,'PV')
    load('PV.dat')
    var = PV;  
  elseif strcmpi(varname,'RH')
    load('RH.dat')
    var = RH * 100;
  elseif strcmpi(varname,'r_eff')
    load('r_eff.dat')
    var = r_eff*1e6;
    load('RH.dat')
    var(RH<1) = 0;
  elseif strcmpi(varname,'LWC')
    load('LWC.dat')
    var = LWC;
    load('RH.dat')
    var(RH<1) = 0;
  elseif strcmpi(varname,'rain_rate')
    load('rain_rate.dat')
    var = rain_rate * 3600*24;
    var(var<0) = 0;
  elseif strcmpi(varname,'S')
    load('RH.dat')
    var = (RH-1)*100;
    var(var<0) = 0;
  elseif strcmpi(varname,'drop_area')
    load('drop_area.dat')
    var = drop_area*1e4;
    var(var<0) = 0;
  elseif strcmpi(varname,'N_drops')
    load('N_drops.dat')
    var = N_drops*1E-6;
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
      diameter(i)=diameter(i-1)*(2.5D-6/diameter(1))^(1D0/(nr_bins-1));
    end
    h = diameter * 1E9;
    size_distrib = true;
  end
 
  ind = time >= 1; 
  ind = 1:length(time);
  [X,Y] = meshgrid(time(ind), h);
  
  
  
  figure;
  surf(X, Y, var(ind,:)','edgecolor','none')
  view(2)
  axis tight
  hc=colorbar;
  cm = colormap;
  set(hc,'fontsize',15)
  set(gca,'color',cm(1,:))
  caxis([quantile(var(:),0.01), quantile(var(:),0.99)])
  %shading interp
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
  varname = regexprep(varname,'_',' ');
  if nargin == 1
    title(varname, 'fontsize', 15)
  else
    title(chemname, 'fontsize', 15)
  end
  set(gca,'fontsize',15)
  if ~size_distrib
    ylim([0, 2000])
  end
  
  
  

end