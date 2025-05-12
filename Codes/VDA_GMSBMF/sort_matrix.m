function res = sort_matrix(score_matrix, interact_matrix)
    % For each disease in columns, sort the score values (predicted)
    % and change the values of the interaction matrix accordingly.
    % In downstream step, go down to each row, compute performance metrics for top 1, 2, 3, .. to summarize the final results.

    sort_index = score_matrix;
    for i = 1:size(sort_index, 2)
        [~, o] = sort(-score_matrix(:, i), 'ascend');
        sort_index(:, i) = o;
    end

    score_sorted = score_matrix;
    y_sorted = interact_matrix;
    for i = 1:size(interact_matrix, 2)
        score_sorted(:, i) = score_matrix(sort_index(:, i), i);
        y_sorted(:, i) = interact_matrix(sort_index(:, i), i);
    end

    res.y_sorted = y_sorted;
    res.score_sorted = score_sorted;
    res.sort_index = sort_index;
end
