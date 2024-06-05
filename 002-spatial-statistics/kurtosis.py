import matplotlib.pyplot as plt
import numpy as np
import scipy.stats as stats

# Set random seed for reproducibility
np.random.seed(42)

# Generate random samples
normal_data = np.random.normal(0, 1, 1000)
t_data = np.random.standard_t(5, 1000)
uniform_data = np.random.uniform(-2, 2, 1000)

# Calculate kurtosis
normal_kurtosis = stats.kurtosis(normal_data, fisher=False)
t_kurtosis = stats.kurtosis(t_data, fisher=False)
uniform_kurtosis = stats.kurtosis(uniform_data, fisher=False)

# Create subplots
fig, axs = plt.subplots(1, 3, figsize=(15, 5))

# Plot Normal Distribution
axs[0].hist(normal_data, bins=30, edgecolor='black', alpha=0.7)
axs[0].set_xlabel('Value')
axs[0].set_ylabel('Frequency')
axs[0].set_title(f'Normal Distribution (Kurtosis: {normal_kurtosis:.2f})')

# Plot T Distribution
axs[1].hist(t_data, bins=30, edgecolor='black', alpha=0.7)
axs[1].set_xlabel('Value')
axs[1].set_ylabel('Frequency')
axs[1].set_title(f'T Distribution (Kurtosis: {t_kurtosis:.2f})')

# Plot Uniform Distribution
axs[2].hist(uniform_data, bins=30, edgecolor='black', alpha=0.7)
axs[2].set_xlabel('Value')
axs[2].set_ylabel('Frequency')
axs[2].set_title(f'Uniform Distribution (Kurtosis: {uniform_kurtosis:.2f})')

# Adjust spacing between subplots
plt.tight_layout()

# Display the plot
plt.show()