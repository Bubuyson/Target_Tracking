function sensor_list = createSensorList(sensor_options, vars)
    fields = fieldnames(sensor_options);
    sensor_list = strings(0);
    for i = 1:length(fields)
        field = fields{i};
        if field == "radar"
            if sensor_options.radar
                sensor_list(end+1) = vars.radar_name;
            end
        else
            str = "Field name " + field + " is not handled";
            error(str);
        end
    end
end