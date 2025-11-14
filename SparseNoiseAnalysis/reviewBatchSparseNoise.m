%% load and concatenate all units

dataDir = 'D:\Code\ephysAnalysis2024\SparseNoiseAnalysis';
dataFiles = dir(fullfile(dataDir, 'M25104_20250902*')); % ''dmTuning_statRun_8t*'dgTuning_statrun

tic
for ifile = 1:numel(dataFiles)
    ifile
    load(fullfile(dataFiles(ifile).folder, dataFiles(ifile).name), 'units');
    % load(fullfile(dataFiles(ifile).folder, dataFiles(ifile).name),  'trials'); % just trials
    session(ifile).units = units;
    % session(ifile).trials = trials;
end
toc

units = cat(1,session.units)

%% subplots 

doplot=true;
if doplot
% idx=1:numel(units);
idx = find([units.firing_rate]>1 & [units.shank]==0);

tic
for iunit = 1:numel(idx)
    % For every new set of panels, open a new figure (except for the first set)
    if mod(iunit-1, 35) == 0
        % if iunit > 1
            % pause;      % Wait for user input after a set of 20
            % close(gcf); % Close the current figure
        % end
        figure;         % Open a new figure
    end

    % Create subplot in a 4x5 grid
    subplot(7, 5, mod(iunit-1, 35) + 1)

    % Plot the receptive field
    thisMap_s = units(idx(iunit)).smthRF;
    imagesc(flipud(thisMap_s))
    xlabel('Azimuth')
    ylabel('Elevation')
    xticks(linspace(1, size(thisMap_s, 2), 13))
    xticklabels(-120:20:120)
    yticks(linspace(1, size(thisMap_s, 1), 7))
    yticklabels(-60:20:60)
    title(sprintf('Unit %i', idx(iunit)))
    colorbar
    set(gca, "TickDir", "out", 'box', 'off', 'Color', 'none', 'FontSize', 8)
    % axis equal,
    axis xy
    abv = abs(max(units(idx(iunit)).smthRF(:)));
    % colormap(redblue), clim([-abv abv])
end

toc
end

%%

good_idx = [40,43,46,53,57,58,63,67,68,73,74,75,80,82,83,87,89,90,91,94,95,96,97,100,102,103,104,106,107,109,110,112,114,115,116,117,120,122,123,124,126,127,131,132,135,136,137,138,140,142,143,144,145,146,148,158,159,160,161,163,165,166,167,169,170,171,172,173,174,175,177,178,181,182,188,193,194,195,196,201,203,205,206,208,209,212,213,215,217,218,220,221,222,223,225,229,235,237,238,239,240,241,242,246,247,249,250,251,252,253,254,255,259,260,261,262,263,264,265,269,270,271,273,274,275,277,278,283,284,285,287,288,296,299,300,301,312,313,314,316,317,324,325,327,329,331,335,341,342,...
            524,525,530,531,532,541,542,544,545,550,555,557,560,561,562,566,577,579,581,585,587,594,595,596,598,599,603,606,607,614,615,617,622,624,625,629,632,633,635,637,638,639,640,641,642,643,647,648,649,650,651,654,655,656,657,659,660,661,662,663,669,673,675,676,678,679,683,684,686,687,692,694,695,697,705,706,708,712,713,714,715,717,718,719,720,725,727,729,738,744,748,750,753,758,760,761,762,764,766,780,895,909,912,915,924,926,932,933,937,940,950,951,952,953,955,956,958,959,961,962,963,964,965,966,967,968,971,973,975,978,979,980,982,...
            699,985,986,1049,1055,1060,1068,1071,1086,1103,1105,1107,1109,1112,1118,1119,1130,1138,1144,1158,1177,1178,1180,1185,1187,1193,1194,1195,1196,1197,1198,1202,1195,1196,1197,1198,1202,1204,1206,1210,1211,1212,1213,1214,1216,1221,1224,1225,1228,1230,1234,1235,1238,1240,1246,1253,1264,1271,1272,1274,1275,1276,1306,1338,1360,1362,1427,1448,1456,1465,1468,1471,1472,1477,1480,1481,1488,1489,1490,1491,1493,1496,1497,1503,1513,1517,1518,1520,1522,1528,1531,1543,1544,1546,1552,1556,1557,1558,...
            1562, 1564 1576 1577 1578 1586 1587 1599 1600 1603 1606 1607 1608 1611 1614 1617 1618 1621 1623 1627 1629 1675 1685 1701 1703 1704 1705 1701 1703 1704 1705 1710 1712 1715 1724 1738 1741 1742 1743 1749 1759 1760 1763 1770 1784 1785 1798 1799 1809 1812 1823 1835 1837 1847 1861 1866 1867 1869 1872 1876 1878 1883 1884 1885 1888 1895 1896 2067 2096 2183 2184 2185 2188 2189 2196 2204 2205 2207 2209 2211 2213 2215 2217 2221 2222 2223 2226 2227 2228 2234 2235 2236 2240 2241 2244 2247 2249 2250 2252 2253 2254 2258 2259 2260 2261 2263 2265 2268 2271 2272 2274 2275 2276 2283 2284 2290 2291 2292 2293 2297]% 4
allocs= cat(1,units.location);

figure, hold on
plot(allocs(:,1), allocs(:,2),'k.');
plot(allocs([good_idx],1), allocs([good_idx],2),'r*')

%%

figure

nDepthSplits=10;
maxVal = max(allocs(:,2));
minVal = min(allocs(:,2));
edges = linspace(minVal,maxVal,nDepthSplits+1);
depth_disc = discretize(allocs(:,2),edges)';

iplot=0;
for idepth = nDepthSplits:-1:1

    for ishank=0:3

        iplot=iplot+1;
        try
            subplot(nDepthSplits,4, iplot)
            idx = find([units.shank]==ishank & depth_disc==idepth);
            rf = nanmean(cat(3,units(idx).normRF),3);
            imagesc(flipud(rf))
            xlabel('Azimuth')
            ylabel('Elevation')
            xticks(linspace(1,size(rf,2),13))
            xticklabels(-120:20:120)
            yticks(linspace(1,size(rf,1),7))
            yticklabels(-60:20:60)
            title(sprintf('Shank %i, Depth %i',ishank,idepth))
            % colorbar
            % colormap(1-gray)
            axis xy
        catch
        end
    end
end