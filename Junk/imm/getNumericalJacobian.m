function H = getNumericalJacobian(fn, pt, my_eps)
    if nargin ~= 2 && nargin ~= 3
        error('Invalid usage of getJacobian')
    end
    if size(pt, 2) ~= 1
        error('Invalid size for getJacobian function')
    end
    if nargin == 3
        eps = my_eps;
    else
        if length(pt) == 3
            eps = [1e-4 1e-5 1e-5]';
        else
            eps = 1e-5 * ones(length(pt), 1);
        end
    end

    if size(eps, 2) ~= 1
        error('Invalid size for getJacobian function')
    end 

    if ~(class(fn) == "function_handle")
        H = fn;
        return
    end

    hx_original = fn(pt);

    H = zeros(length(hx_original), length(pt));
    for i = 1:length(pt)
        perturbed_pt = pt;
        perturbed_pt(i) = perturbed_pt(i) + eps(i);
        hx_perturbed = fn(perturbed_pt);
        H(:, i) = (hx_perturbed - hx_original) / eps(i);
    end
    % perturbed_pts = repmat(pt, 1, length(hx_original)) + diag(eps);
    % hx_perturbed = fn(perturbed_pts);
    % 
    % H = (hx_perturbed - hx_original) ./ my_eps;


end