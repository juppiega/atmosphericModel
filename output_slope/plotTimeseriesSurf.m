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
    name = 'r_{eff} [\mu m]';
    var = onlyCloud(var, RH);
  elseif strcmpi(varname,'LWC')
    load('LWC.dat')
    var = LWC*1e3;
    name = 'LWC [g/m^3]';
    load('RH.dat')
    var = onlyCloud(var, RH);
  elseif strcmpi(varname,'rain_rate')
    load('rain_rate.dat')
    var = rain_rate * 3600*24;
    var(var<0) = 0;
    name = 'Rain rate [mm/day]';
  elseif strcmpi(varname,'S')
    load('RH.dat')
    var = (RH-1)*100;
    var(var<0) = 0;
  elseif strcmpi(varname,'drop_area')
    load('drop_area.dat')
    var = drop_area*1e4;
    var(var<0) = 0;
    name = 'Drop area [cm^2/m^3]';
  elseif strcmpi(varname,'N_drops')
    load('N_drops.dat')
    var = N_drops*1E-6;
    name = 'Drop concentration [1/cm^3]';
  elseif strcmpi(varname,'theta')
    load('theta.dat')
    var = theta;
    name = '\Theta [K]';
  elseif strcmpi(varname,'chemistry')
    var = load([chemname,'.dat']);
    if strcmpi(chemname, 'size_distribution')
      var = log10(var/1E6);
      var(var < 0) = 0;
    end  
    var = var(time >= 2, :);
    time = time(time >= 2);
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
    name = 'Kh [m^2/s]';
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
 
  ind = find(time >= 2); 
  %ind = 1:length(time);
  [X,Y] = meshgrid((0:time(2)-time(1):1.001)*24, h);
  
  N_day = round(1 / (time(2)-time(1))) + 1;
  %v = cat(3, var(ind(1)+(0:N_day-1),:), var(ind(1)+N_day+(0:N_day-1),:));
  %var = mean(v,3);
  
  if strcmpi(varname,'w_kin')
      var(var < 0) = 0;
      name = 'w_u [m/s]';
  end
  
  if strcmpi(varname,'flux_t')
      name = '\Theta flux [K m/s]';
  end
  
  figure;
  set(gcf,'color','w');
  surf(X, Y, var(ind(1):ind(end)-1,:)','edgecolor','none')
  view(2)
  axis tight
  colorbar
  caxis([quantile(var(:),0.01), quantile(var(:),0.99)])
  %shading interp
  if  strcmpi(varname,'flux_t')
      v = var(ind(1):ind(end)-1,:);
      caxis([min(v(:)), -min(v(:))])
  end
  if size_distrib
    set(gca, 'YScale', 'log')
  end
  xlabel('Time [h]', 'fontsize', 15)
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
  if ~exist('name', 'var')
      title(varname, 'fontsize', 15)
  else
      title(name,'fontsize',15)
  end
  set(gca,'fontsize',15)
  if ~size_distrib
    ylim([0, 3000])
  end
  
  
  

end

function var = onlyCloud(var, RH)

[LCL, top] = find_LCL_top(RH);

for i = 1:length(LCL)
    v = var(i,LCL(i):top(i));
    var(i,:) = 0;
    var(i,LCL(i):top(i)) = v;
end

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