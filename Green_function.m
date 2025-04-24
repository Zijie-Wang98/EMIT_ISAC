function G = Green_function(k,p_vec)
    p = vecnorm(p_vec);
    p_dir = p_vec/p;
    g = exp(1j*k*p)/(4*pi*p);
    G = g*((1-1j/(k*p)-1/(k^2*p^2))*eye(3)-(1-3*1j/(k*p)-3/(k^2*p^2))*(p_dir*p_dir.'));
    G = G(:,1:2);
end