figure()
set(gcf, 'Units', 'point', 'Position', [200 200 1000, 300], 'color', 'w')
fontsize = 12;
strengths_to_plot = 1:length(alphas);


%----------------------------------------
% sequenceness vs simulation alpha power
%----------------------------------------
subplot(1,3,1)
[~, id] = min(abs(fft_freq_domain - 10));
alpha_psd = squeeze(mean(psd_store(:,strengths_to_plot,:,id),3)); % [exp, alpha];
shadedErrorBar(alphas(strengths_to_plot), mean(alpha_psd), std(alpha_psd)./sqrt(size(alpha_psd,1)), 'k', 0);
xlabel({'MEG alpha power (a.u.)', ''})
ylabel('Sequenceness alpha power (PSD)')
title({'MEG \alpha power and sequenceness periodicity', ''})
box off
set(gca, 'FontSize', fontsize)


%----------------------------------------
% ground truth sequenceness at 50ms, normalized effect (as a proportion of the group perm threshold)
%----------------------------------------
subplot(1,3,2)
pk = 5;
grp_peak_effect = squeeze(mean(sqn_store(:,strengths_to_plot,:,pk), 3));        % [exp, alpha]
norm_grp_pk = grp_peak_effect./GroupThresh_store(:, strengths_to_plot);         % [exp, alpha]
shadedErrorBar(alphas(strengths_to_plot), mean(norm_grp_pk), std(norm_grp_pk)/sqrt(size(norm_grp_pk,1))); hold on
plot([alphas(strengths_to_plot(1)) alphas(strengths_to_plot(end))], [1 1], 'r:', 'LineWidth', 5)
xlabel({'MEG alpha power (a.u.)', ''})
ylabel('Sequenceness effect (normalized)')
title({'MEG \alpha power and sequenceness SNR', ''})
box off
set(gca, 'FontSize',  fontsize)


%----------------------------------------
% PSD by alpha
%----------------------------------------
subplot(1,3,3)
h = [];
frange = 3:15;
psd_by_alpha = squeeze(mean(psd_store(:,strengths_to_plot,:,frange),3)); % [exp, alpha, freq]
h = [];
for aa = 1:length(strengths_to_plot)
    h(aa) = plot(fft_freq_domain(frange), mean(squeeze(psd_by_alpha(:,aa,:))), 'LineWidth', 3, 'Color', [45*(aa) 0 200]/256);
    hold on
end

hh = {};
for ii = 1:length(h)
    hh{ii} = sprintf('\\alpha = %d', alphas(strengths_to_plot(ii)));
end
legend(h, hh)
ylabel('PSD (mean over sim.)')
xlabel('Frequency (Hz)')
title({'Sequenceness PSD by MEG alpha', ''})
box off
set(gca, 'FontSize', fontsize)