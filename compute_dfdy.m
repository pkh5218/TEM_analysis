function dfdy = compute_dfdy(f, phi)
% f: [Ny x Np], phi: latitude (rad)
    dphi = phi(2)-phi(1);
    dfdy = zeros(size(f));
    dfdy(2:end-1,:) = (f(3:end,:) - f(1:end-2,:)) / (2*dphi);
    dfdy(1,:)       = (f(2,:) - f(1,:)) / dphi;
    dfdy(end,:)     = (f(end,:) - f(end-1,:)) / dphi;
end