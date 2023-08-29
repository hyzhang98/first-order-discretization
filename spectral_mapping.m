function [Y, obj, obj2] = spectral_mapping(F, L, D, eta)
    % 
    % F (n * c): Relaxed solution
    % L (n * n): Laplacian matrix
    % D (n * n): Degree matrix; it should be eye(n) if using Ratio Cut
    % eta: recommend 10**-4 by default 
    % Output:
    % Y (n * c): Discrete solution
    % obj: ||F * R - G||^2 + eta * <Delta, 2LFR>
    % obj2: tr(G' * L * G)
    % Author: Hongyuan Zhang
    
    [n, n_clusters] = size(F);
    max_iter = 10;
    % Initialize R
    R = eye(n_clusters);
    % Initialize Y
    Y = rand(n, n_clusters);
    y = max(Y, [], 2);
    Y = double(Y == y);
    obj = zeros(max_iter, 1);
    obj2 = zeros(max_iter, 1);
    for itr = 1:max_iter
        % Solve Y
        M = F * R - eta * L * F * R;
        d_tmp = diag(sqrt(D)) .* M;
        for i = 1:n 
            d = sum(d_tmp .* Y, 1);
            % ns = sum(Y.^2, 1);
            ns = diag((Y' .* diag(D)') * Y)';

            inc = d + sqrt(D(i, i)) * M(i, :) .* (1 - Y(i, :));
            inc = inc ./ sqrt(ns + D(i,i) * (1 - Y(i, :))); 

            base = d - sqrt(D(i, i)) * M(i, :) .* Y(i, :);
            base = base ./ sqrt(ns - D(i,i) * Y(i, :)); 

            gain = inc - base; 
            [~, idx] = max(gain);
            new_yi = zeros(1, n_clusters);
            new_yi(idx) = 1;
            old_yi = Y(i, :);
            Y(i, :) = new_yi; 
            [~, y] = max(Y, [], 2);
            if length(unique(y)) ~= n_clusters
                Y(i, :) = old_yi;
            end
        end

        % Solve R
        % diag(Y' * Y);
        % T = (Y' * D * Y ).^0.5;
        d = diag(D);
        T = ((Y .* d)' * Y).^0.5;
        t = diag(T);
        t(t~=0) = t(t~=0).^-1;
        % G = D^0.5 * Y * diag(t); 
        G = Y .* d.^0.5 * diag(t);
        T = F' * G - eta * F' * L * G;
        [U, ~, V] = svd(T);
        oldR = R;
        R = U * V';
        if norm(oldR - R, 'fro') < 10^-3
            %break;
        end
        % obj(itr) = norm(F * R - G, 'fro')^2 + 2 * eta * trace(R' * F' * L * (G-F*R));
        % obj2(itr) = trace(G' * L * G);
    end
    
end