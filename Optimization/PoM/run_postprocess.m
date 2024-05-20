% enter the name of the folder you want to postprocess
%typically of the format outputyyyymmddTxxxxx/

filename = 'output20220910T150241/';

%Enter the number of runs you want to analyze.
%E.g. if you had 1000 cells enter 1000 for rawsize
rawsize = 50;

%For the calculation of Activity fraction and deltapeaks
init_volt = -63;
nbins = 100;



logConstraint = zeros(rawsize,7);
count = 0;

for i = 1:rawsize
    file = strcat(filename,num2str(i));
    load(file);
    AP = real(Y(:,14));
    
    %Number of peaks 
    pks = findpeaks(AP);

    if numel(pks) >= 100
        count = count + 1;
        %Activity Fraction
        activityFraction = (numel(find(AP>init_volt)))/numel(AP);
        %Delta Peak
        [APcounts,edges] = histcounts(AP,nbins);
        [max1, ind1] = max(APcounts);
        APcounts(ind1) = -Inf;
        [max2, ind2] = max(APcounts);
        deltapeak = (max1 - max2)/max1;
    
        %Peaks based on FFT
        APfft = nufft(AP,T);
        pks1 = findpeaks(abs(APfft));

        %Mean Cytosolic Calcium
        cytCal = mean(real(Y(:,13)));
        peakcytCal = max(real(Y(:,13)));

        %Peaks of Cytosolic Calcium based on FFT
        Cafft = nufft(real(Y(:,13)),T);
        pks2 = findpeaks(abs(Cafft));
        logConstraint(count,1) = activityFraction;
        logConstraint(count,2) = deltapeak;
        logConstraint(count,3) = numel(pks);
        logConstraint(count,4) = numel(pks1);
        logConstraint(count,5) = cytCal;
        logConstraint(count,6) = peakcytCal;
        logConstraint(count,7) = numel(pks2);
    end
end

if count < rawsize
    logConstraint(count+1:rawsize,:) = [];
end

%Find busters
idx = find(((logConstraint(:,4)>20000)&(logConstraint(:,4) <= 40000))~=0);

if length(idx) ~=0
    bursters = logConstraint(idx,:);
end

%Find spikers
idx = find(((logConstraint(:,4)>40000))~=0);
if length(idx) ~=0
    spikers = logConstraint(idx,:);
end

%Find plateau
idx = find(((logConstraint(:,4) <= 20000))~=0);
if length(idx) ~=0
    plateau = logConstraint(idx,:);
end
    
%write a xlsx file with:
%column 1: activity fraction
%column 2: deltapeak
%column 3: Number of peaks in the voltage signal
%column 4: Number of peaks in the FFT of the voltage
%column 5: Mean of cytosolic calcium
%column 6: Peak of cytosolic calcium
%column 7: Number of peaks in the FFT of the calcium
writematrix(spikers,'Classes.xlsx','Sheet',1);
writematrix(bursters,'Classes.xlsx','Sheet',2);
writematrix(plateau,'Classes.xlsx','Sheet',3);

