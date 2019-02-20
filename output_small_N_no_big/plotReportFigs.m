function [] = plotReportFigs ()

fs = 21;

q_k = load('q_k.dat');
q_m = load('q_m.dat');
flux_q_k = load('flux_q_k.dat');
flux_q_m = load('flux_q_m.dat');
flux_local = load('flux_q_local_m.dat');
flux_mass = flux_q_m - flux_local;
%flux_mass(abs(flux_mass) < 1E-6) = 0;

time = load('time.dat');
load('h.dat');

q_rel = q_k./q_m - 1;
flux_rel = flux_q_k./flux_q_m - 1;

ind = find(time >= 3.5); 
time = time * 24;
[X,Y] = meshgrid(time(ind) - time(ind(1)) + 12, h);

figure;
subplot(2,3,1)
surf(X, Y, q_m(ind,:)','edgecolor','none')
view(2)
axis tight
colorbar('fontsize',fs)
xlim([12, 24])
c_m = caxis;
xlabel('time [h]','fontsize',fs);
ylabel('z [m]','fontsize',fs);
title('q_{NL}','fontsize',fs)
set(gca,'fontsize',fs)
set(gca, 'color', 'b')

subplot(2,3,2)
surf(X, Y, q_k(ind,:)','edgecolor','none')
view(2)
axis tight
colorbar('fontsize',fs)
xlim([12, 24])
caxis(c_m);
xlabel('time [h]','fontsize',fs);
ylabel('z [m]','fontsize',fs);
title('q_L','fontsize',fs)
set(gca,'fontsize',fs)
set(gca, 'color', 'b')

subplot(2,3,3)
surf(X, Y, q_rel(ind,:)','edgecolor','none')
view(2)
axis tight
colorbar('fontsize',fs)
caxis([-.1,.1])
xlim([12, 24])
xlabel('time [h]','fontsize',fs);
ylabel('z [m]','fontsize',fs);
title('q_L / q_{NL} - 1','fontsize',fs)
set(gca,'fontsize',fs)
set(gca, 'color', 'b')

[X,Y] = meshgrid(time(ind) - time(ind(1)) + 12, 0.5*(h(1:end-1) + h(2:end)));

subplot(2,3,4)
surf(X, Y, flux_q_m(ind-1,:)','edgecolor','none')
view(2)
axis tight
colorbar('fontsize',fs)
xlim([12, 24])
c_k = caxis;
xlabel('time [h]','fontsize',fs);
ylabel('z [m]','fontsize',fs);
title('Flux q_{NL}','fontsize',fs)
set(gca,'fontsize',fs)
set(gca, 'color', 'b')

subplot(2,3,5)
surf(X, Y, flux_q_k(ind-1,:)','edgecolor','none')
view(2)
axis tight
colorbar('fontsize',fs)
xlim([12, 24])
caxis(c_k);
xlabel('time [h]','fontsize',fs);
ylabel('z [m]','fontsize',fs);
title('Flux q_L','fontsize',fs)
set(gca,'fontsize',fs)
set(gca, 'color', 'b')

subplot(2,3,6)
flux_mass_rel = max(min(flux_mass ./ (flux_mass + flux_local),1),0);
Z = flux_rel(ind-1,:);
Z(flux_mass_rel(ind-1,:) == 0 & flux_q_k(ind-1,:) < 1E-5) = 0;
surf(X, Y, Z','edgecolor','none')
view(2)
axis tight
colorbar('fontsize',fs)
caxis([-1,1])
xlim([12, 24])
xlabel('time [h]','fontsize',fs);
ylabel('z [m]','fontsize',fs);
title('Flux q_L / Flux q_{NL} - 1','fontsize',fs)
set(gca,'fontsize',fs)
set(gca, 'color', 'b')

start_ind = ind(1);
%[X,Y] = meshgrid(time(ind) - time(ind(1)) + 12, 0.5*(h(1:end-1) + h(2:end)));
%figure;
%surf(X, Y, flux_rel(start_ind-1:end,:)','edgecolor','none')
%view(2)
%axis tight
%colorbar('fontsize',fs)
%caxis([-2.5E-4,2.5E-4])
%xlim([12, 24])
%title('flux_K - flux_{NL}','fontsize',fs)
%xlabel('time [h]','fontsize',fs);
%ylabel('z [m]','fontsize',fs);
%set(gca,'fontsize',fs)

figure;
t = 17;
[~,closest_ind] = min(abs(time/24 - (3+t/24)));
plot(q_k(closest_ind,:),h,'linewidth',2.0,q_m(closest_ind,:),h,'linewidth',2.0)
ylabel('z [m]','fontsize',fs);
legend('q_L','q_{NL}')
set(gca,'fontsize',fs)
title('Local (L), nonlocal (NL) q at 1800h','fontsize',fs)

figure;
i = start_ind-1:size(flux_rel,1);
flux_mass_rel = max(min(flux_mass ./ (flux_mass + flux_local),1),0);
surf(X, Y, flux_mass_rel(i,:)','edgecolor','none')
view(2)
axis tight
colorbar('fontsize',fs)
%caxis([0,1])
xlim([12, 24])
title('Mass flux to total flux','fontsize',fs)
xlabel('time [h]','fontsize',fs);
ylabel('z [m]','fontsize',fs);
set(gca,'fontsize',fs)
set (gca, 'color', 'b')

%[X,Y] = meshgrid(time(ind) - time(ind(1)) + 12, h(2:end-1));
%
%figure;
%i = start_ind-1:size(flux_rel,1);
%dz = h(2)-h(1);
%fm_conv = -diff(flux_mass,1,2)/dz;
%fl_conv = -diff(flux_local,1,2)/dz;
%flux_mass_rel = fm_conv ./ (fl_conv);
%surf(X, Y, flux_mass_rel(i,:)','edgecolor','none')
%view(2)
%axis tight
%colorbar('fontsize',fs)
%%caxis([0,1])
%xlim([12, 24])
%title('Mass flux convg. to total flux convg.','fontsize',fs)
%xlabel('time [h]','fontsize',fs);
%ylabel('z [m]','fontsize',fs);
%set(gca,'fontsize',fs)
%set (gca, 'color', 'b')

%figure;
%t = 17;
%[~,closest_ind] = min(abs(time/24 - (3+t/24)));
%hmid = 0.5*(h(1:end-1) + h(2:end));
%load('hd.dat');
%y = flux_mass_rel(closest_ind,:); y(hmid > hd(closest_ind)) = 0;
%plot(y, hmid,'linewidth',2.0);
%title('Mass flux to total flux of q','fontsize',fs)
%xlabel('Ratio','fontsize',fs);
%ylabel('z [m]','fontsize',fs);
%set(gca,'fontsize',fs)



[X,Y] = meshgrid(time(ind) - time(ind(1)) + 12, h);

figure;
subplot(3,2,1)
Z = load('theta.dat');
surf(X, Y, Z(ind(1):end,:)','edgecolor','none')
view(2)
axis tight
colorbar('fontsize',fs)
caxis([293,297])
title('θ [K]','fontsize',fs)
%xlabel('time [h]','fontsize',fs);
ylabel('z [m]','fontsize',fs);
set(gca,'fontsize',fs)
xlim([12, 24])
set (gca, 'color', 'b')

subplot(3,2,5)
Z = load('mass_flux.dat');
surf(X, Y, Z(ind(1)-1:end,:)','edgecolor','none')
view(2)
axis tight
colorbar('fontsize',fs)
title('Mass flux [m/s]','fontsize',fs)
xlabel('Time [h]','fontsize',fs);
ylabel('z [m]','fontsize',fs);
set(gca,'fontsize',fs)
xlim([12, 24])
set (gca, 'color', 'b')

[X,Y] = meshgrid(time(ind) - time(ind(1)) + 12, 0.5*(h(1:end-1) + h(2:end)));

subplot(3,2,2)
Z = load('flux_t.dat');
surf(X, Y, Z(ind(1)-1:end,:)','edgecolor','none')
view(2)
axis tight
colorbar('fontsize',fs)
caxis([min(Z(:)),-min(Z(:))])
title('θ flux [K m/s]','fontsize',fs)
%xlabel('time [h]','fontsize',fs);
%ylabel('z [m]','fontsize',fs);
set(gca,'fontsize',fs)
xlim([12, 24])
set (gca, 'color', 'b')

subplot(3,2,3)
Z = load('E_tot.dat');
surf(X, Y, Z(ind(1)-1:end,:)','edgecolor','none')
view(2)
axis tight
colorbar('fontsize',fs)
title('E_{tot} [m^2 / s^2]','fontsize',fs)
%xlabel('time [h]','fontsize',fs);
ylabel('z [m]','fontsize',fs);
set(gca,'fontsize',fs)
xlim([12, 24])
set (gca, 'color', 'b')

subplot(3,2,4)
Z = load('Kh.dat');
surf(X, Y, Z(ind(1)-1:end,:)','edgecolor','none')
view(2)
axis tight
colorbar('fontsize',fs)
%caxis([0,200])
title('K_h [m^2 / s]','fontsize',fs)
%xlabel('time [h]','fontsize',fs);
%ylabel('z [m]','fontsize',fs);
set(gca,'fontsize',fs)
xlim([12, 24])
set (gca, 'color', 'b')

subplot(3,2,6)
Z = load('mixing_length.dat');
surf(X, Y, Z(ind(1)-1:end,:)','edgecolor','none')
view(2)
axis tight
colorbar('fontsize',fs)
title('Mixing length [m]','fontsize',fs)
xlabel('Time [h]','fontsize',fs);
%ylabel('z [m]','fontsize',fs);
set(gca,'fontsize',fs)
xlim([12, 24])
set (gca, 'color', 'b')


end
