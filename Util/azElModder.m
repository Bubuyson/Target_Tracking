function new_az_el = azElModder(az_el, option)
    if nargin < 2
        mult = 180 / pi;
        inv_mult = pi/180;
    else
        if option == "degree"
            mult = 1;
            inv_mult = 1;
        else
            mult = 180 / pi;
            inv_mult = pi/180;
        end
    end
    az = az_el(1, :) * mult;
    el = az_el(2, :) * mult;

    new_az_el = zeros(2, size(az_el, 2));
    for i = 1:size(az_el, 2)
        if az(i) > 180
            if el > 90
                az(i) = azModder(az(i));
                el(i) = elPositiveModder(el(i));
            elseif el(i) < -90
                az(i) = azModder(az(i));
                el(i) = elNegativeModder(el(i));
            else
                az(i) = azModder(az(i));
            end
        elseif az(i) < -180
            if el(i) > 90
                az(i) = azModder(az(i));
                el(i) = elPositiveModder(el(i));
            elseif el(i) < -90
                az(i) = azModder(az(i));
                el(i) = elNegativeModder(el(i));
            else
                az(i) = azModder(az(i));
            end
        else
            if el(i) > 90
                if az(i) > 0
                    az(i) = az(i) - 180;
                else
                    az(i) = az(i) + 180;
                end
                el(i) = elPositiveModder(el(i));
            elseif el(i) < -90
                if az(i) > 0
                    az(i) = az(i) - 180;
                else
                    az(i) = az(i) + 180;
                end
                el(i) = elNegativeModder(el(i));
            else
            end
        end
        new_az_el(:, i) = [az(i) * inv_mult;
                           el(i) * inv_mult];
    end

    

    function az = azModder(az)
        az = mod(az + 180, 360) - 180;
    end

    function el = elPositiveModder(el)
        el = -mod(el, 90) + 90;
    end

    function el = elNegativeModder(el)
        el = -mod(el, 90);
    end

    % % % function makeUnitTest()                                             %#ok
    % % %     test_values = [150 91; 150 -91; 181 60; 181 91; -181 60; -181 -91; 60 60; 181 -91; -181 91]';
    % % %     disp('Testing test value: ')
    % % %     disp(test_values)
    % % %     val = azElModder(test_value);
    % % %     disp('Output value: ')
    % % %     disp(val)
    % % % end
end