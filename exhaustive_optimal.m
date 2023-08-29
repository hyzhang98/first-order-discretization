function [y_pred, obj] = exhaustive_optimal(L, D, n_clusters)
    n = size(L, 1); 
    Y = zeros(n, n_clusters);
    [obj_left, Y_left] = exhaustive_Y(Y, L, D, 1, 0);
    [obj_right, Y_right] = exhaustive_Y(Y, L, D, 1, 1);
    [obj, idx] = min([obj_left, obj_right]);
    if idx == 1
        Y_pred = Y_left;
    else
        Y_pred = Y_right;
    end
    [~, y_pred] = max(Y_pred, [], 2);
end

function [obj, Y_pred] = exhaustive_Y(Y, L, D, i, left)
    Y(i, 1) = left;
    Y(i, 2) = 1-left;
    n = size(L, 1);
    if n == i 
        if sum(Y(:, 1)) == n || sum(Y(:, 2)) == n
            obj = 100;
            Y_pred = Y;
            return;
        end
        obj = cal_obj(Y, L, D);
        Y_pred = Y;
        return
    end
    [obj_left, Y_left] = exhaustive_Y(Y, L, D, i+1, 1);
    [obj_right, Y_right] = exhaustive_Y(Y, L, D, i+1, 0);
    [obj, idx] = min([obj_left, obj_right]);
    if idx == 1
        Y_pred = Y_left;
    else
        Y_pred = Y_right;
    end
end