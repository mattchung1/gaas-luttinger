% Helper function for building the Hamiltonian

function H = buildHamiltonian(kx, ky, kz, p)
    % p: struct of constants (hbar, m0, gamma1, gamma2, gamma3, Delta)

    Q = (p.hbar^2 * p.gamma2 / (2*p.m0)) * (kx^2 + ky^2 - 2*kz^2);
    R = (p.hbar^2 / (2*p.m0)) * (-sqrt(3)*p.gamma2*(kx^2 - ky^2) + 1i*2*sqrt(3)*p.gamma3*kx*ky);
    S = (p.hbar^2 * p.gamma3 / p.m0) * sqrt(3) * (kx - 1i*ky) * kz;
    P = (p.hbar^2 * p.gamma1 / (2*p.m0)) * (kx^2 + ky^2 + kz^2);

    H = zeros(6,6);
    H(1,1) = -P - Q;
    H(2,2) = -P + Q;
    H(3,3) = -P + Q;
    H(4,4) = -P - Q;
    H(5,5) = -P - p.Delta;
    H(6,6) = -P - p.Delta;

    H(1,2) = S;
    H(1,3) = -R;
    H(1,5) = S/sqrt(2);
    H(1,6) = -sqrt(2)*R;

    H(2,1) = conj(S);
    H(2,4) = -R;
    H(2,5) = sqrt(2)*Q;
    H(2,6) = -sqrt(1.5)*S;

    H(3,1) = -conj(R);
    H(3,4) = -S;
    H(3,5) = -sqrt(1.5)*conj(S);
    H(3,6) = -sqrt(2)*Q;

    H(4,2) = -conj(R);
    H(4,3) = -conj(S);
    H(4,5) = sqrt(2)*conj(R);
    H(4,6) = conj(S)/sqrt(2);

    H(5,1) = conj(S)/sqrt(2);
    H(5,2) = sqrt(2)*conj(Q);
    H(5,3) = -sqrt(1.5)*S;
    H(5,4) = sqrt(2)*R;

    H(6,1) = -sqrt(2)*conj(R);
    H(6,2) = -sqrt(1.5)*conj(S);
    H(6,3) = -sqrt(2)*conj(Q);
    H(6,4) = S/sqrt(2);
end