function innov = twoPiMod(y, y_, index)
    innov = y - y_;
    modded_val = mod(innov(index) + pi, 2* pi) - pi;
    innov(index) = modded_val;
end