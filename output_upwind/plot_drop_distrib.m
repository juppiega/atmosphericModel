function plot_drop_distrib()

upwind = load('drop_distribution.upwind.dat'); 
slope = load('drop_distribution.slope.dat'); 

ind = 6:2:size(slope,1);
time = 0.25*(1:size(slope,1)) - 1; 

N_bins = length(upwind);
max_diameter = 3250D-6;
diameter(1)=4D-6;
for i=2:N_bins
    diameter(i)=diameter(i-1)*(max_diameter/diameter(1))^(1D0/(N_bins-1));
end
midpoints = zeros(size(diameter));
midpoints(2:end) = exp(0.5*(log(diameter(1:end-1)) + log(diameter(2:end))));
midpoints(1) = midpoints(2) - (midpoints(3) - midpoints(2));

figure;
for i = 1:length(ind)
    subplot(2,2,i)
    n_up = upwind(ind(i),1:end-1)*1e-3./ (midpoints(2:end) - midpoints(1:end-1));
    n_sl = slope(ind(i),1:end-1)*1e-3./ (midpoints(2:end) - midpoints(1:end-1));
    n_sl(diameter(1:end-1)>1E-3) = 0;
    loglog(diameter(1:end-1)*1e3, n_up, diameter(1:end-1)*1e3, n_sl, 'linewidth', 2.0)
    xlabel('D_p [mm]', 'fontsize', 15)
    ylabel('N(D_p) [m^{-3} mm^{-1}]','fontsize', 15)
    legend('Upwind', 'PLM')
    set(gca,'fontsize',15)
    title(['t: ',num2str(time(ind(i))),' h'],'fontsize',15)
    ylim([1e-10,max(n_sl(:))])
    xlim([diameter(1), diameter(end-1)]*1e3)
end

end