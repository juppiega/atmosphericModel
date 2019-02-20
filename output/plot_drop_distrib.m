function plot_drop_distrib()

upwind = load('drop_distribution.upwind.dat'); upwind = upwind(end,:);
slope = load('drop_distribution.slope.dat'); slope = slope(end,:);

N_bins = length(upwind);
max_diameter = 3250D-6;
diameter(1)=4D-6;
for i=2:N_bins
    diameter(i)=diameter(i-1)*(max_diameter/diameter(1))^(1D0/(N_bins-1));
end
midpoints = zeros(size(diameter));
midpoints(2:end) = exp(0.5*(log(diameter(1:end-1)) + log(diameter(2:end))));
midpoints(1) = midpoints(2) - (midpoints(3) - midpoints(2));
n_up = upwind(1:end-1)*1e-3./ (midpoints(2:end) - midpoints(1:end-1));
n_sl = slope(1:end-1)*1e-3./ (midpoints(2:end) - midpoints(1:end-1));

figure;
loglog(diameter(1:end-1)*1e3, n_up, diameter(1:end-1)*1e3, n_sl, 'linewidth', 2.0)
xlabel('D_p [mm]', 'fontsize', 15)
ylabel('N(D_p) [m^{-3} mm^{-1}]','fontsize', 15)
legend('Upwind', 'LPM')
set(gca,'fontsize',15)
title('Drop distribution','fontsize',15)
ylim([1e-10,max(n_sl(:))])

end