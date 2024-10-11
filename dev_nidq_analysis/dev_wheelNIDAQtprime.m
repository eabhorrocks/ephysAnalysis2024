%% get wheel data

wheel_a_up = table2array(readtable(fullfile(nidqDir,'wheel_a_up_tprime.txt')));
wheel_a_down = table2array(readtable(fullfile(nidqDir,'wheel_a_down_tprime.txt')));
wheel_b_up = table2array(readtable(fullfile(nidqDir,'wheel_b_up_tprime.txt')));
wheel_b_down = table2array(readtable(fullfile(nidqDir,'wheel_b_down_tprime.txt')));


maxVal = max([wheel_a_up; wheel_a_down; wheel_b_up; wheel_b_down]);

% get ttl pulses that actually correspond to tick changes
wheel_a_valid=wheel_a_up;
tic
for ival = 1:numel(wheel_a_up)
    val = wheel_a_up(ival);
    temp=val-wheel_b_up;
    temp=temp(temp>0);
    [bup, iup]=min(temp);
    temp=val-wheel_b_down;
    temp=temp(temp>0);
    [bdown, idown]=min(temp);
    if bdown>bup
        wheel_a_valid(ival) = nan; % remove if not valid
    end
end
time2finisha = toc;

% get ttl pulses that actually correspond to tick changes
wheel_b_valid=wheel_b_up;
tic
for ival = 1:numel(wheel_b_up)
    val = wheel_b_up(ival);
    temp=val-wheel_a_up;
    temp=temp(temp>0);
    [aup, iup]=min(temp);
    temp=val-wheel_a_down;
    temp=temp(temp>0);
    [adown, idown]=min(temp);
    if adown>aup
        wheel_b_valid(ival) = nan;
    end
end
time2finishb = toc;

%% remove invalid values

wheel_a_valid(isnan(wheel_a_valid))=[];
wheel_b_valid(isnan(wheel_b_valid))=[];

%% resample the tick position

resampleRate=100;

num_a = numel(wheel_a_valid);
num_b = numel(wheel_b_valid);

[tickChangeTimes, idx] = sort(cat(1,wheel_a_valid(:),wheel_b_valid(:)));

tickChanges = ones(numel(idx),1);
tickChanges(idx>num_a)=-1;
tickPos = cumsum(tickChanges);

resampleTime = 0:(1/resampleRate):ceil(max(tickChangeTimes));
resampleTickPos = interp1(tickChangeTimes,tickPos,resampleTime);
