function [class, outs, param] = Andrean(pks,troughs,V,T,Ca,param)
    outs = cell(6,1);
    if isempty(pks)
        class = 'Silent';
        outs{1} = NaN;
        outs{2} = NaN;
        outs{3} = NaN;
        outs{4} = NaN;
        outs{5} = mean(Ca);
        outs{6} = max(Ca);
    else 
        if length(troughs)>length(pks)
            average = mean(V(pks)-V(troughs(1:length(pks))));
        elseif length(troughs)<length(pks)
            if length(pks)-length(troughs) > 1
                average = mean(V(pks(1:length(troughs)))-V(troughs));
            else
                average = mean(V(pks(1:end-1))-V(troughs));
            end  
        else 
            average = mean(V(pks)-V(troughs));
        end
        if average > 10
            [crosses,~] = peakfinder(V, 10, -50, -1, 1);
            if isempty(crosses)
                class = 'Depolarized';
                outs{1} = NaN;
                outs{2} = NaN;
                outs{3} = NaN;
                outs{4} = NaN;
                outs{5} = mean(Ca);
                outs{6} = max(Ca);
            elseif length(pks) < 2
                class = 'Depolarized';
                outs{1} = NaN;
                outs{2} = NaN;
                outs{3} = NaN;
                outs{4} = NaN;
                outs{5} = mean(Ca);
                outs{6} = max(Ca);
            else
                maxint = max(diff(T(pks)));
                if maxint > 2000 && numel(pks) > 15 % previously 25 when simulations were 1000s long (now 600s)
                    diffs = diff(T(pks));
                    end_pks = pks(find(diffs>2000));
                    start_pks = zeros(length(end_pks)+1,1);
                    end_pks(end+1) = pks(end);
                    for i = 2:length(end_pks)
                        start_pks(i) = pks(find(pks==end_pks(i-1))+1);
                    end
                    start_pks(1) = pks(1);
                    for i = 1:length(start_pks)
                        burst_dur(i) = mean(T(end_pks(i))-T(start_pks(i)));
                        startind = find(pks==start_pks(i));
                        endind = find(pks==end_pks(i));
                        burst_freq(i) = 1000/mean(diff(T(pks(startind:endind))));
                    end
                    num_bursts = length(start_pks);
                    burst_dur = mean(burst_dur);
                    burst_freq = mean(burst_freq);                
                    class = 'Bursting';
                    outs{1} = burst_freq;
                    outs{2} = burst_dur;
                    outs{3} = num_bursts;
                    outs{4} = NaN;
                    outs{5} = mean(Ca);
                    outs{6} = max(Ca);
                else
                    if numel(pks) < 30 % previously 50 when simulations were 1000s long (now 600s)
                        class = 'Other';
                        outs{1} = NaN;
                        outs{2} = NaN;
                        outs{3} = NaN;
                        outs{4} = NaN;
                        outs{5} = mean(Ca);
                        outs{6} = max(Ca);
                    else
                        spike_freq = 1000/mean(diff(T(pks)));
                        class = 'Spiking';
                        outs{1} = NaN;
                        outs{2} = NaN;
                        outs{3} = NaN;
                        outs{4} = spike_freq;
                        outs{5} = mean(Ca);
                        outs{6} = max(Ca);
                    end
                end
            end
        end
    end
end