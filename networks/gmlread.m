function A = gmlread(fileName)

fid = fopen(fileName);
stringAsc = fread(fid)';
fclose(fid);
string = char(stringAsc);
isDirect = ~isempty(regexp(string,'directed 1\n','once'));

patternNode = 'node[^a-z].*?]';
[strStart,strStop] = regexp(string,patternNode);
stringTmp = string(strStart(1):strStop(1));
if isempty(regexp(stringTmp,'label', 'once'))
    isLabel = 0;
    nodeList(1:length(strStart)) = struct('id',[]);
else
    isLabel = 1;
    nodeList(1:length(strStart)) = struct('id',[],'label',[]);
end
for i = 1 : length(strStart)
    stringTmp = string(strStart(i):strStop(i));
    [strStartTmp,strStopTmp] = regexp(stringTmp,'id [0-9]+\n');
    nodeList(i).id = str2double(stringTmp(strStartTmp+3:strStopTmp-1));
    if isLabel
        [strStartTmp,strStopTmp] = regexp(stringTmp,'label ".*"\n');
        nodeList(i).label = stringTmp(strStartTmp+7:strStopTmp-2);
    end
end

patternEdge = 'edge[^a-z].*?]';
[strStart,strStop] = regexp(string,patternEdge);
edgeList = zeros(length(strStart),2);
for i = 1 : length(strStart)
    stringTmp = string(strStart(i):strStop(i));
    [strStartTmp,strStopTmp] = regexp(stringTmp,'source [0-9]+\n');
    edgeList(i,1) = str2double(stringTmp(strStartTmp+7:strStopTmp-1));
    [strStartTmp,strStopTmp] = regexp(stringTmp,'target [0-9]+\n');
    edgeList(i,2) = str2double(stringTmp(strStartTmp+7:strStopTmp-1));
end

A = zeros(length(nodeList));
for i = 1 : length(edgeList)
    A(edgeList(i,1)+1-nodeList(1).id,edgeList(i,2)+1-nodeList(1).id) = 1;
    if isDirect == 0
        A(edgeList(i,2)+1-nodeList(1).id,edgeList(i,1)+1-nodeList(1).id) = 1;
    end
end
end