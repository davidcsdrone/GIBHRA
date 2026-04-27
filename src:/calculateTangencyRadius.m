function[r_new, closest_pair_ids] = calculateTangencyRadius(triplet_ids, x_loc, y_loc)
%Calculates the uniform radius required to make the closest pair of nodes
%in a given triplet mutually tangent.

    id_i = triplet_ids(1);
    id_j = triplet_ids(2);
    id_k = triplet_ids(3);

    C_i = [x_loc(id_i), y_loc(id_i)];
    C_j = [x_loc(id_j), y_loc(id_j)];
    C_k = [x_loc(id_k), y_loc(id_k)];

    d_ij = sqrt(sum((C_i - C_j).^2));
    d_ik = sqrt(sum((C_i - C_k).^2));
    d_jk = sqrt(sum((C_j - C_k).^2));


    distances = [d_ij, d_ik, d_jk];
    [d_min, min_idx] = min(distances);

    if min_idx == 1
        closest_pair_ids = [id_i, id_j];
    elseif min_idx == 2
        closest_pair_ids = [id_i, id_k];
    else
        closest_pair_ids = [id_j, id_k];
    end
    r_new = d_min / 2;
end