function declare_info_level(scope, level)
% InfoLevel is declared globally, so we can use it without check within any class that
% have this declaration.
%
% useage: 
%   (1) declare default information level (|scope| and |level| is not declared). If the
%       corresponding information level has been assigned a value, the default value will
%       not take effect. 
%   (2) set information level with |scope| and |level|. This will override existing value.
global InfoLevel;
if nargin >= 2
    if ischar(scope)
        scope = {scope};
    end
    for i = 1:length(scope)
        InfoLevel.(scope{i}) = level(i);
    end
else
    default_level = struct('Global', DisplayLevel.Off, ...
        'Class', DisplayLevel.Final, ...
        'ClassDebug', DisplayLevel.Notify, ...
        'InnerModel', DisplayLevel.Off, ...
        'InnerModelDebug', DisplayLevel.Notify, ...
        'UserModel', DisplayLevel.Final, ...
        'UserModelDebug', DisplayLevel.Notify);
    InfoLevel = structmerge(InfoLevel, default_level, 'exclude');
end
end

