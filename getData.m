function output = getData (pathName)
rawXData = xlsread(pathName);
rawXData = transpose(rawXData);
xData = [];
for k = 1:length(rawXData)
    xData = [xData; rawXData(k)];
end
output = xData;
end

