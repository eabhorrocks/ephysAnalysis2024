%% get set of intervals for trials
statTrials = allFlickerTrials([allFlickerTrials.runFlag]==0);
runTrials = allFlickerTrials([allFlickerTrials.runFlag]==1);

for ifreq = 1:numel(flickFreqs)
    % stationary/locomotion trials
    stimTrials = statTrials([statTrials.StimulusFrequency]==flickFreqs(ifreq));
    stat_intervalsCell{ifreq}= [vertcat(stimTrials.start_time), vertcat(stimTrials.stop_time)]*1000;

    stimTrials = runTrials([runTrials.StimulusFrequency]==flickFreqs(ifreq));
    run_intervalsCell{ifreq}= [vertcat(stimTrials.start_time), vertcat(stimTrials.stop_time)]*1000;
end



%% chronux settings

params.Fs=1000;
params.fpass = [0 150];
params.trialave = 0;
params.pad = 0;
params.tapers = [2 4];
params.err =[1 0.05];

%%
for iunit = 1%:numel(units)
    close
    figure
for ifreq = 1:12;
try
    %spikeTimes = units(iunit).spike_times*1000;
    spikeTimes = cat(1,units.spike_times)*1000;

    tints = stat_intervalsCell{ifreq};
    trialLength = round(mean(tints(:,2)-tints(:,1)),-1,'decimals'); % round to 10ms

    stimOnTimes = tints(:,1);
    stimOffTimes = tints(:,2);
    intStarts = stimOnTimes;
    intStops = stimOnTimes+trialLength;

    ipsth = 1;
    % loop through intervals (usually trials) and bin spikes
    clear p stat_t

    for itrial = 1:size(tints,1)

        % spikes that occur during this interval
        relevantSpikeTimes = spikeTimes(spikeTimes>=intStarts(itrial) &...
            spikeTimes <= intStops(itrial));
        % 0-centre to stimulus onset. save these for later raster plot
        p(ipsth).trial(itrial).trialCentricSpikes = relevantSpikeTimes - stimOnTimes(itrial);
        p(ipsth).trial(itrial).nSpikes = numel(p(ipsth).trial(itrial).trialCentricSpikes);
        stat_t(itrial).spike_times = p(1).trial(itrial).trialCentricSpikes/1000;

    end






    [S,f,R,Serr]=mtspectrumpt(stat_t,params);
    stat_t(itrial).S = S; stat_t(itrial).f=f; stat_t(itrial).R=R; stat_t(itrial).Serr = Serr;
    stat_S=S; stat_f = f; stat_R = R; stat_Serr = Serr;
    stat_S = pow2db(stat_S);

    % run

    % tints = run_intervalsCell{ifreq};
    % trialLength = round(mean(tints(:,2)-tints(:,1)),-1,'decimals'); % round to 10ms
    % 
    % stimOnTimes = tints(:,1);
    % stimOffTimes = tints(:,2);
    % intStarts = stimOnTimes;
    % intStops = stimOnTimes+trialLength;
    % 
    % ipsth = 1;
    % % loop through intervals (usually trials) and bin spikes
    % clear p
    % 
    % for itrial = 1:size(tints,1)
    % 
    %     % spikes that occur during this interval
    %     relevantSpikeTimes = spikeTimes(spikeTimes>=intStarts(itrial) &...
    %         spikeTimes <= intStops(itrial));
    %     % 0-centre to stimulus onset. save these for later raster plot
    %     p(ipsth).trial(itrial).trialCentricSpikes = relevantSpikeTimes - stimOnTimes(itrial);
    %     p(ipsth).trial(itrial).nSpikes = numel(p(ipsth).trial(itrial).trialCentricSpikes);
    %     run_t(itrial).spike_times = p(1).trial(itrial).trialCentricSpikes/1000;

    % end





    % 
    % [S,f,R,Serr]=mtspectrumpt(run_t,params);
    % 
    % run_t(itrial).S = S; run_t(itrial).f=f; run_t(itrial).R=R; run_t(itrial).Serr = Serr;
    % run_S=S; run_f = f; run_R = R; run_Serr = Serr;
    % run_S = pow2db(run_S);



    subplot(2,6,ifreq)
    shadedErrorBar(f,mean(stat_S,2),sem(stat_S,2))
    % shadedErrorBar(f,mean(run_S,2),sem(run_S,2),'lineprops','r')
    ax=gca; ax.XTick = [1, 10:10:150]; grid on
    ylim([25 50])
    title(flickFreqs(ifreq))
catch
end
end


iunit

pause


end


