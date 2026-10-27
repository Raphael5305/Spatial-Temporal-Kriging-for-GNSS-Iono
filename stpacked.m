function out = stpacked(packed, s, t)
    hasLUT = isfield(packed,'s_grid') && ~isempty(packed.s_grid) && ...
             isfield(packed,'t_grid') && ~isempty(packed.t_grid) && ...
             isfield(packed,'S_on_s') && ~isempty(packed.S_on_s) && ...
             isfield(packed,'T_on_t') && ~isempty(packed.T_on_t);
    if hasLUT
        out = eval_lut(packed, s, t);
    else
        out = eval_matern2(packed, s, t);
    end
end

function gamma_st = eval_lut(pk, s, t)
    ys = interp1(double(pk.s_grid), double(pk.S_on_s), double(s), 'linear','extrap');
    yt = interp1(double(pk.t_grid), double(pk.T_on_t), double(t), 'linear','extrap');
    gamma_st = ys + yt - double(pk.alpha) .* ys .* yt;
end

function gamma_st = eval_matern2(pk, s, t)
    ys = matern_inline2(double(s), pk.Sparam);
    yt = matern_inline2(double(t), pk.Tparam);
    gamma_st = ys + yt - double(pk.alpha) .* ys .* yt;
end

function g = matern_inline2(d, P)
    d   = max(d, 0);
    rho = max(double(P.rho), eps);
    nu  = double(P.nu);
    s2  = double(P.partialsill);
    c0  = double(P.nugget);
    g   = zeros(size(d), 'like', d);
    mask = (d > 0);
    if ~any(mask), return; end
    r = d(mask) ./ rho;
    if abs(nu-0.5) < 1e-12
        rho_mt = exp(-r);
    elseif abs(nu-1.5) < 1e-12
        c = sqrt(3);  rho_mt = (1 + c*r) .* exp(-c*r);
    elseif abs(nu-2.5) < 1e-12
        c = sqrt(5);  rho_mt = (1 + c*r + (5/3)*(r.^2)) .* exp(-c*r);
    else
        coef   = 2^(1-nu) / gamma(nu);
        rho_mt = coef .* (r.^nu) .* besselk(nu, r);
    end
    g(mask) = c0 + s2 .* (1 - rho_mt);
end
