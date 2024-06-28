import pandas as pd
import numpy as np
from scipy import stats
from sklearn.model_selection import LeaveOneOut
from sklearn.metrics import accuracy_score
from astropy.io import fits
from astropy.table import Table

def save_to_fits(data, filename):
    astropy_table = Table.from_pandas(data)
    astropy_table.write(filename, format='fits', overwrite=True)


def load_data(stars_path, quasars_path):
    stars = pd.read_csv(stars_path)
    quasars = pd.read_csv(quasars_path)
    return stars, quasars

def remove_contaminated(stars, quasars):
    # Assuming 'objID' in stars and 'specObjID' in quasars are supposed to be equivalent identifiers
    # Remove any star entries that have an objID appearing in quasars' specObjID
    clean_stars = stars[~stars['objID'].isin(quasars['specObjId'])]
    return clean_stars

def calculate_colors(df):
    df['u-g'] = df['psfmag_u'] - df['psfmag_g']
    df['g-r'] = df['psfmag_g'] - df['psfmag_r']
    df['r-i'] = df['psfmag_r'] - df['psfmag_i']
    df['i-z'] = df['psfmag_i'] - df['psfmag_z']
    return df[['u-g', 'g-r', 'r-i', 'i-z']]

def fit_kde(data, bandwidth):
    # Transpose the data to fit KDE: expects (n_features, n_samples)
    values = data.T
    kde = stats.gaussian_kde(values, bw_method=bandwidth)
    return kde

def calculate_posterior(kde_star, kde_quasar, data_point, prior_star, prior_quasar):
    prob_star = np.exp(kde_star.logpdf(data_point))
    prob_quasar = np.exp(kde_quasar.logpdf(data_point))
    
    posterior_star = (prob_star * prior_star) / (prob_star * prior_star + prob_quasar * prior_quasar)
    
    return posterior_star < 0.5

# def cross_validate_bandwidth(data, labels, bandwidths):
#     loo = LeaveOneOut()
#     best_accuracy = 0
#     best_bandwidth = None
    
#     for bandwidth in bandwidths:
#         accuracies = []
#         for train_index, test_index in loo.split(data):
#             train_data, test_data = data.iloc[train_index], data.iloc[test_index]
#             train_labels, test_labels = labels[train_index], labels[test_index]
#             kde = fit_kde(train_data, bandwidth)
#             predictions = [calculate_posterior(kde, kde, test_data.iloc[i], 0.88, 0.12) for i in range(len(test_data))]
#             accuracies.append(accuracy_score(test_labels, predictions))
        
#         mean_accuracy = np.mean(accuracies)
#         if mean_accuracy > best_accuracy:
#             best_accuracy = mean_accuracy
#             best_bandwidth = bandwidth
    
#     return best_bandwidth, best_accuracy

def main():
    stars_path = '/hpc/group/cosmology/denniswu/scolnic/stars.csv'
    quasars_path = '/hpc/group/cosmology/denniswu/scolnic/quasars.csv'
    
    stars, quasars = load_data(stars_path, quasars_path)

    # Remove contaminated quasars from the stars catalog
    clean_stars = remove_contaminated(stars, quasars)

    print("finished remove")

    stars_colors = calculate_colors(clean_stars)
    quasars_colors = calculate_colors(quasars)

    print("finished colors")
    
    # Leave-one-out cross-validation for bandwidth selection
    bandwidths = np.linspace(0.1, 0.2, 5)
    best_bandwidth_star = 0.1 # cross_validate_bandwidth(stars_colors[:10000], np.zeros(len(stars_colors[:10000])), bandwidths)
    best_bandwidth_quasar = 0.15 # cross_validate_bandwidth(quasars_colors[:1000], np.ones(len(quasars_colors[:1000])), bandwidths)

    print("finished leave one out")
    
    # Fit KDE models with the best bandwidths found
    kde_star = fit_kde(stars_colors, best_bandwidth_star)
    kde_quasar = fit_kde(quasars_colors, best_bandwidth_quasar)

    print("finished fitting")
    
    # Vectorized calculation of posteriors
    data_points = np.vstack([stars_colors[col].values for col in ['u-g', 'g-r', 'r-i', 'i-z']])  # Prepare data for vectorized computation
    prob_star = np.exp(kde_star.logpdf(data_points))  # Probability of each point under the star KDE
    prob_quasar = np.exp(kde_quasar.logpdf(data_points))  # Probability of each point under the quasar KDE
    posterior_star = prob_star * 0.88 / (prob_star * 0.88 + prob_quasar * 0.12)  # Posterior probability under star class
    is_quasar = posterior_star < 0.5  # Quasar classification
    stars['is_quasar'] = is_quasar

    # Cut in redshift space
    new_quasars = stars[(stars['is_quasar']) & (stars['z'] > 1) & (stars['z'] < 2.2)]
    print("finished classify")
    
    # Assuming the DataFrame stars includes columns 'ra', 'dec', and magnitudes
    required_columns = ['ra', 'dec', 'psfmag_u', 'psfmag_g', 'psfmag_r', 'psfmag_i', 'psfmag_z']
    new_quasars = new_quasars[required_columns]
    
    # Save new quasars to a FITS file
    save_to_fits(new_quasars, 'new_quasars.fits')

if __name__ == "__main__":
    main()
