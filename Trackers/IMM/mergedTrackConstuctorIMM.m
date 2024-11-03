function track = mergedTrackConstuctorIMM(varargin)
    if nargin == 0
        track = struct('track_id', {}, 'x', {}, 'P', {});
    else
        field_names = varargin(1:2:end);
        field_values = varargin(2:2:end);
        if length(field_values) ~= length(field_names)
            field_names(end) = [];
        end
        track = struct();
        for i = 1:length(field_names)
            track.(field_names{i}) = field_values{i};
        end
    end
    
end