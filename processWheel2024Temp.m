function wheel = processWheel2024Temp(wheelFile,smoothMethod,windowSize)

wheel_temp = readtable(wheelFile);
wheel.pos = wheel_temp.Wheel;
wheel.time = wheel_temp.ArduinoTime;
halfMax = max(wheel.pos)/2;
wheel_unwrapped = unwrap(wheel.pos, halfMax);
wheel.dist = wheel2unit(wheel_unwrapped, 1024, 19);%.*-1; % pos, ticks/rev, wheel diam

wheel.rawSpeed = diff(wheel.dist)./diff(wheel.time);
wheel.rawSpeed = movmean(wheel.rawSpeed, 2);
wheel.rawSpeed = [wheel.rawSpeed(1); wheel.rawSpeed];

windowSize_bin = windowSize/10; %10ms arduino sample rate
wheel.smthSpeed = smoothdata(wheel.rawSpeed,smoothMethod,windowSize_bin);
wheel.smthSpeed = wheel.smthSpeed * 1000; % convert from cm/ms to cm/s


end