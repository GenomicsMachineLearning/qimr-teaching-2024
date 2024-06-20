import numpy as np
import matplotlib.pyplot as plt


def moran_i(values, weights):
    n = len(values)
    mean_value = np.mean(values)
    norm = np.sum(weights)

    numerator = 0
    denominator = 0
    for i in range(n):
        for j in range(n):
            numerator += weights[i][j] * (values[i] - mean_value) * (values[j] - mean_value)
        denominator += (values[i] - mean_value) ** 2

    moran = (n / norm) * (numerator / denominator)
    return moran


# Create a sample dataset
values = [1, 2, 3, 4, 5]
coordinates = [(1.2, 1.0), (1.5, 1.5), (2.2, 2.6), (3.0, 3.3), (3.5, 3.5)]
#coordinates = [(1.2, 1.0), (1.5, 1.1), (2.2, 2.6), (3.0, 3.3), (3.1, 3.5)]
#coordinates = [(3.1, 3.5), (1.2, 1.0), (1.5, 1.1), (3.0, 3.3), (2.2, 2.6)]

# Calculate the squared inverse distance weights
weights = np.zeros((len(values), len(values)))
for i in range(len(values)):
    for j in range(len(values)):
        if i != j:
            dist = np.sqrt((coordinates[i][0] - coordinates[j][0]) ** 2 + (coordinates[i][1] - coordinates[j][1]) ** 2)
            weights[i][j] = 1 / (dist ** 2)

# Calculate Moran's I
moran = moran_i(values, weights)

# Visualize the data points and their relationships
plt.figure(figsize=(6, 6))
plt.scatter([coord[0] for coord in coordinates], [coord[1] for coord in coordinates], c=values, cmap='viridis', s=100)

for i in range(len(values)):
    for j in range(len(values)):
        if weights[i][j] > 0:
            plt.plot([coordinates[i][0], coordinates[j][0]], [coordinates[i][1], coordinates[j][1]], 'k-',
                     linewidth=0.5)
            mid_x = (coordinates[i][0] + coordinates[j][0]) / 2
            mid_y = (coordinates[i][1] + coordinates[j][1]) / 2
            dist = np.sqrt((coordinates[i][0] - coordinates[j][0]) ** 2 + (coordinates[i][1] - coordinates[j][1]) ** 2)
            plt.text(mid_x, mid_y, f"{1 / dist ** 2:.2f}", fontsize=12, ha='center', va='center')

plt.colorbar(label='Value')
plt.xlabel('X')
plt.ylabel('Y')
plt.title(f"Moran's I: {moran:.3f}")
plt.show()