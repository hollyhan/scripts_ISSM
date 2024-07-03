function perturbed_ice_thickness = perturb_ice_thickness_data(ice_thickness)
    % This function takes in a matrix of ice thickness data and returns a perturbed version of it
    % ice_thickness: Input ice thickness data (number of grid points x number of time steps)

    % Set the random number generator for reproducibility
    rng(123); %comment out if want to produce more than one random perturbed field

    % Analyze the existing ice thickness field
    disp('     Getting covariance of the thickness field')
    [mean_thickness, empirical_covariance] = analyze_ice_thickness(ice_thickness);
    
    % Generate Gaussian Random Field (GRF) with the same covariance structure
    n_grid_points = size(ice_thickness, 1);
    n_time_steps = size(ice_thickness, 2);
    disp('     Generating Gaussian Random Field')
    grf = generate_grf_with_covariance(n_grid_points, n_time_steps, empirical_covariance);
    
    % Perturb the existing ice thickness field
    perturbed_ice_thickness = ice_thickness + grf;

    % Function to analyze the existing ice thickness field
    function [mean_thickness, empirical_covariance] = analyze_ice_thickness(ice_thickness)
        % Calculate the mean thickness for each grid point across time steps
        mean_thickness = mean(ice_thickness, 2);
        
        % Center the ice thickness field by subtracting the mean for each grid point
        centered_thickness = ice_thickness - mean_thickness;
        
        % Compute the empirical covariance matrix across grid points
        empirical_covariance = cov(centered_thickness');
        
        % Regularize the covariance matrix
        epsilon = 1e-9; % Small regularization term. Make it bigger if fails (1e-10 was too small)
        empirical_covariance = empirical_covariance + epsilon * eye(size(empirical_covariance));
    end
    
    % Function to generate a Gaussian Random Field (GRF) for each time step
    function grf = generate_grf_with_covariance(n_grid_points, n_time_steps, covariance_matrix)
        % Generate white noise for each grid point and each time step
        white_noise = randn(n_grid_points, n_time_steps);  % Generate noise for each grid point, for each time step
        
        % Generate GRF by multiplying white noise with the Cholesky decomposition of the covariance matrix
        chol_covariance = chol(covariance_matrix, 'lower');
        grf = chol_covariance * white_noise;
    end
end


