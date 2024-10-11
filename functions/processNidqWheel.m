function wheel = processNidqWheel(wheel,samplingFreq,smoothMethod,windowSize)

% resample wheel signal specified frequency
sampleSpacing = 1/samplingFreq;
xq = 0:sampleSpacing:ceil(max(wheel.time));
wheel.pos = interp1(wheel.time, wheel.pos,xq); 
wheel.pos(1)=0; 
wheel.time = xq; 

% make column vectors
wheel.pos=wheel.pos(:); wheel.time=wheel.time(:);

% convert from ticks to distance in metric units (normally cm)
wheel.dist = wheel2unit(wheel.pos, 1024, 19);%.*-1; % pos, ticks/rev, wheel diam

% calculate the speed by taking differential w/ respect to time
wheel.rawSpeed = diff(wheel.dist)./diff(wheel.time);
wheel.rawSpeed = movmean(wheel.rawSpeed, 2);
wheel.rawSpeed = [wheel.rawSpeed(1); wheel.rawSpeed];

% find the smoothing window in binsizes and smooth raw speed signal
windowSize_bin = windowSize/(1*sampleSpacing);
wheel.smthSpeed = smoothdata(wheel.rawSpeed,smoothMethod,windowSize_bin);

end