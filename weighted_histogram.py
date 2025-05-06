import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# === Filepath to Excel file ===
filepath = "C:/Users/wharrick25/OneDrive - Charlotte Country Day School/torsionGLUFullSet.xlsx"

# === Read Excel file into a DataFrame ===
df = pd.read_excel(filepath)

# === Convert DataFrame to 2D matrix ===
data_matrix = df.to_numpy()

def calculate_integral(df_hist, bins, lower_bound, upper_bound):
    """
    Calculate the integral of the histogram between two points.
    """
    # Get bin heights and edges
    hist, bin_edges = np.histogram(df_hist["value"], bins=bins, weights=df_hist["occupancy"])
    
    # Calculate the integral
    integral = 0
    for i in range(len(bin_edges) - 1):
        if bin_edges[i] >= lower_bound and bin_edges[i + 1] <= upper_bound:
            integral += hist[i] * (bin_edges[i + 1] - bin_edges[i])
    return integral

for columnno in range(3, 16, 2):
    column = []
    bondname = data_matrix[0][columnno]
    for i in range(2, len(data_matrix) - 1):
        column.append((data_matrix[i][columnno], data_matrix[i][columnno + 1]))

    # === Create DataFrame for histogram ===
    df_hist = pd.DataFrame(column, columns=["value", "occupancy"])

    # === Ensure numeric data ===
    df_hist["value"] = pd.to_numeric(df_hist["value"], errors="coerce")
    df_hist["occupancy"] = pd.to_numeric(df_hist["occupancy"], errors="coerce")
    df_hist = df_hist.dropna()  # Drop rows with NaN values

    # === Calculate bin size and edges dynamically ===
    numberOfBins = 100  # Restore the original number of bins
    data_min, data_max = df_hist["value"].min(), df_hist["value"].max()
    bin_size = (data_max - data_min) / numberOfBins  # Adjust the divisor for desired granularity
    bins = np.arange(data_min, data_max + bin_size, bin_size)

    # === Create weighted histogram ===
    plt.hist(df_hist["value"], bins=bins, weights=df_hist["occupancy"], edgecolor="black")

    # === Add gridlines for each bin ===
    plt.grid(True, which="both", axis="x", linestyle="--", linewidth=0.5)  # Add gridlines for all bins

    # === Set x-axis ticks and labels ===
    plt.xticks(bins[::5], rotation=45)  # Label every 10th bin

    # === Label plot ===
    plt.title(bondname)
    plt.xlabel("Angle")
    plt.ylabel("# of Occurrences")

    # === Show/save plot ===
    plt.tight_layout()
    bondname = bondname.replace(" ", "")  
    plt.savefig(f"histoimages/{bondname}histo.png", dpi=300)
    plt.clf()  # Clear the current figure for the next plot

    # === Calculate and print the integral for a specific range ===
    #if bondname=="psi":
    #    lower_bound = 140  # Replace with your desired lower bound
    #    upper_bound = 220  # Replace with your desired upper bound
    #    integral = calculate_integral(df_hist, bins, lower_bound, upper_bound)
    #    print(f"Integral of {bondname} between {lower_bound} and {upper_bound}: {integral}")
