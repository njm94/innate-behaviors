function numComponents = pca_image(image, N)
    % Input:
    % image: HxW video where H is height, W is width
    %
    % Output:
    % numComponents: The number of principal components that explain N% of the variance

%     image_reshaped = image(:);
% 
%     % If video contains nans, remove those
%     image_reshaped(isnan(image_reshaped)) = [];

    % Perform PCA (using 'econ' to calculate only the necessary number of components)
    [~, ~, latent] = pca(image, 'economy', true); % Transpose for PCA

    % Compute the cumulative explained variance
    explained_variance = cumsum(latent) / sum(latent);

    % Find the number of components that explain 90% of the variance
    numComponents = find(explained_variance >= N/100, 1);
end