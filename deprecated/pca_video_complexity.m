function numComponents = pca_video_complexity(video, N)
    % Input:
    % video: HxWxT video where H is height, W is width, T is number of frames
    % N:     Percentage of variance trying to achieve
    %
    % Output:
    % numComponents: The number of principal components that explain 90% of the variance

    % Get video dimensions
    [H, W, T] = size(video);

    % Reshape the video to 2D: each column is a flattened frame
    video_reshaped = reshape(video, H*W, T);

    % Perform PCA (using 'econ' to calculate only the necessary number of components)
    [~, ~, latent] = pca(video_reshaped', 'economy', true); % Transpose for PCA

    % Compute the cumulative explained variance
    explained_variance = cumsum(latent) / sum(latent);

    % Find the number of components that explain 90% of the variance
    numComponents = find(explained_variance >= N/100, 1);
end