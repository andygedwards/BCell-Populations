function [msg] = plotTS(lowT,highT,V_low,V_high,pks_low,pks_high,trs_low,trs_high,cell_index,plot_title,plot_peaks_and_troughs)
   try
        F = figure; set(F, 'Position', [100, 100, 1200, 1000]);
        subplot(2,1,1), plot(lowT/1000,V_low), hold on, 
        if plot_peaks_and_troughs == 1
            plot(lowT(pks_low)/1000, V_low(pks_low),'r*', lowT(trs_low)/1000, V_low(trs_low),'g*')
        end
        ylim([-100,30])
        ylabel('V_{m} (mV)')
        title('2 mM glucose', 'Fontsize',22)
        subplot(2,1,2), plot(highT/1000,V_high), hold on, 
        if plot_peaks_and_troughs == 1
            plot(highT(pks_high)/1000, V_high(pks_high),'r*', highT(trs_high)/1000, V_high(trs_high),'g*')
        end
        ylim([-100,30])
        ylabel('V_{m} (mV)')
        xlabel('time (s)')
        title('20 mM glucose', 'Fontsize',22)
        sgtitle(F,['Cell number: ', num2str(cell_index),'. ', plot_title], 'Fontsize',24)
        pause
        msg = 'success';
        close all
   catch
        msg = ['Failed while attempting to plot ', 'Cell number ', num2str(cell_index),':', title];
   end
end