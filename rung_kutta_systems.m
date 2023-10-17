function [t, w] = rung_kutta_systems (f, a, b, N, alp)

    % set h, t, and the initial values 
    h = (b - a) / N;
    t = a:h:b;
    w(:, 1) = alp;

    for i = 2 : N + 1
        K1 = h * f(t(i - 1), w(:, i - 1));
        K2 = h * f(t(i - 1) + h / 2, w(:, i - 1) + K1 / 2);
        K3 = h * f(t(i - 1) + h / 2, w(:, i - 1) + K2 / 2);
        K4 = h * f(t(i - 1) + h, w(:, i - 1) + K3);  

        w(:, i) = w(:, i - 1) + (K1 + 2 * K2 + 2 * K3 + K4) / 6;
    end
end