function [stationary_intervals, locomotion_intervals] = ...
    changepointsWheelAnalysis(wheel_zscoreSpeed, wheelSpeed, wheelTime, zThresh, sampleRate,...
    smoothWin, timeBetween, minDuration, minLocoSpeed, maxStatSpeed,...
    preBuffer, postBuffer)

% enforce column vectors
wheel_zscoreSpeed = wheel_zscoreSpeed(:);
wheelSpeed = wheelSpeed(:);
wheelTime = wheelTime(:);


% from Lohani et al., 
[h1,OnTStamp ,OffTStamp ] =changepoints(wheel_zscoreSpeed, zThresh,wheelTime,sampleRate,...
    smoothWin,0,timeBetween,false,false,false);

%% get intitial putative locomotion/stationary epochs

locomotion_intervals = [OnTStamp, OffTStamp];

% use to get putative stationary
end_points = locomotion_intervals(1:end-1, 2);
start_points = locomotion_intervals(2:end, 1);

stationary_intervals = [end_points, start_points];
stationary_intervals = [0, locomotion_intervals(1,1); stationary_intervals]; % add a putative initial interval
stationary_intervals = [stationary_intervals; locomotion_intervals(end,2), max(wheelTime)];

%% remove intervals that don't meet criteria

min_duration=minDuration; % 5s
% min_locoSpeed = 3;
% max_statSpeed = 1;

locoInterval_meanSpeed=[];
locoInterval_duration=[];
locoInterval_timeSinceLast(1)=nan;

for i=1:size(locomotion_intervals,1)
    locoInterval_meanSpeed(i) = mean(wheelSpeed(find(wheelTime==locomotion_intervals(i,1)):find(wheelTime==locomotion_intervals(i,2))));
    locoInterval_duration(i) = locomotion_intervals(i,2)-locomotion_intervals(i,1);
    if i>1
        locoInterval_timeSinceLast(i) = locomotion_intervals(i,1)-locomotion_intervals(i-1,2);
    end
end

loco_to_del = find([locoInterval_meanSpeed]<minLocoSpeed | [locoInterval_duration]<minDuration);

statInterval_meanSpeed=[];
statInterval_duration=[];
statInterval_timeSinceLast(1)=nan;

for i=1:size(stationary_intervals,1)
    statInterval_meanSpeed(i) = mean(wheelSpeed(find(wheelTime==stationary_intervals(i,1)):find(wheelTime==stationary_intervals(i,2))));
    statInterval_duration(i) = stationary_intervals(i,2)-stationary_intervals(i,1);
    if i>1
        statInterval_timeSinceLast(i) = stationary_intervals(i,1)-stationary_intervals(i-1,2);
    end
end


stat_to_del = find([statInterval_meanSpeed]>maxStatSpeed | [statInterval_duration]<min_duration);


locomotion_intervals_filtered = locomotion_intervals;
locomotion_intervals_filtered(loco_to_del,:)=[];

stationary_intervals_filtered = stationary_intervals;
stationary_intervals_filtered(stat_to_del,:)=[];

% remove a buffer period from the start and end of each interval
% pre_buffer=0.5;
% post_buffer=0.5;

stationary_intervals_filtered(:,1) = stationary_intervals_filtered(:,1)+preBuffer;
stationary_intervals_filtered(:,2) = stationary_intervals_filtered(:,2)-postBuffer;

locomotion_intervals_filtered(:,1) = locomotion_intervals_filtered(:,1)+preBuffer;
locomotion_intervals_filtered(:,2) = locomotion_intervals_filtered(:,2)-postBuffer;

stationary_intervals=stationary_intervals_filtered;
locomotion_intervals=locomotion_intervals_filtered;

end