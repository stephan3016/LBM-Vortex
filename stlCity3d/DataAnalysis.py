import numpy as np
import pandas as pd
from matplotlib import pyplot as plt

# Data integration and selection
excel_file_path = r"BaselWind10And100m.xlsx"       #Excel file containing the data set
sheet_name = "Basel"
useful_columns = [1, 3]
skip_rows = range(9)
df = pd.read_excel(excel_file_path, sheet_name=sheet_name, usecols=useful_columns, skiprows=skip_rows)  # Loading useful rows and colums
df=df.dropna()      # Data cleaning -->Getting rid of empty rows


# Reconstructing the distribution:

# First convert everything to m/s instead of km/h
df.iloc[:,0] = df.iloc[:, 0]*1/3600*1000
df.iloc[:,1] = df.iloc[:,1]*1/3600*1000

# The assumption is that both dimensions are independent


# Fitting a logarithmic distribution based on this data
from scipy.stats import lognorm

# Calculate the histogram
hist_10m, bin_edges_10m = np.histogram(df.iloc[:, 0], bins=200)
hist_100m, bin_edges_100m = np.histogram(df.iloc[:, 1], bins=200)

# Calculate the bin centers
bin_centers_10m = (bin_edges_10m[:-1] + bin_edges_10m[1:]) / 2
bin_centers_100m = (bin_edges_100m[:-1] + bin_edges_100m[1:]) / 2

# Calculate the total number of data points
total_count_10m = len(df.iloc[:, 0])
total_count_100m = len(df.iloc[:, 1])

# Calculate the relative frequencies
relative_freq_10m = hist_10m / total_count_10m
relative_freq_100m = hist_100m / total_count_100m

# Fit a logarithmic distribution to the data
params_10m = lognorm.fit(df.iloc[:, 0])
params_100m = lognorm.fit(df.iloc[:, 1])
print(params_100m)

# Generate the probability density function (PDF) values based on the fitted logarithmic distribution
pdf_10m = lognorm.pdf(bin_centers_10m, *params_10m)
pdf_100m = lognorm.pdf(bin_centers_100m, *params_100m)

# Plot the histogram and the fitted logarithmic distribution
fig2, axes2 = plt.subplots()
axes2.hist(df.iloc[:, 0], bins=200, density=True, alpha=0.5, label="z = 100 m", weights=np.ones_like(df.iloc[:, 0]) / total_count_10m)
axes2.plot(bin_centers_10m, pdf_10m, 'b-', label='Logarithmic Fit')
axes2.set_xlabel("Wind velocity magnitude [m/s]")
axes2.set_ylabel("Relative Frequency [-]")

axes2.hist(df.iloc[:, 1], bins=200, density=True, alpha=0.5, label="z = 100 m", weights=np.ones_like(df.iloc[:, 1]) / total_count_100m)
axes2.plot(bin_centers_100m, pdf_100m, 'r-', label='Logarithmic Fit')


# Making both plots comparable:
axes2.set_xlim(min(min(df.iloc[:,0]),min(df.iloc[:,1])), max(max(df.iloc[:,0]),max(df.iloc[:,1])))
axes2.legend()
plt.show()



# Fitting formula: v(z) = a*(z/z_max)^(1/alpha)
# For now just take the expected value
# More fancy statistics can be applied though

#selecting points
Velexpected_10m = lognorm.mean(*params_10m)
Velexpected_100m = lognorm.mean(*params_100m)

#regression (very basic for now, must be improved)
from scipy.optimize import curve_fit

def function(x,a,alpha):
    return a*(x/120)**(1/alpha)

x = [10, 100]
y = [Velexpected_10m, Velexpected_100m]

coefficients, coeffficients_covariance = curve_fit(function, x, y)
print(coefficients)
a_fit, alpha_fit = coefficients


fig3, axes3 = plt.subplots()
X = np.linspace(0,120,120)
Y = []
for i in range(0, len(X)):
    Y.append(function(X[i], a_fit, alpha_fit))

axes3.scatter(y, x)
axes3.plot(Y, X, "b-", label="Fit function")
axes3.set_ylabel("Height [m]")
axes3.set_xlabel("Velocity magnitude [m/s]")
plt.legend()
plt.show()
