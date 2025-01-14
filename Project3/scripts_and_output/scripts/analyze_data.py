import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import os

# Load the data
file_path = "Project3/scripts/final_results.csv"
data = pd.read_csv(file_path)

# Data Cleaning
data = data[data["status"] == "OK"]  # Filter out failed instances
data["randomization"] = data["randomization"].astype(bool)

##########################
# 1) Summaries by Category
##########################
category_summary = data.groupby("category").agg({
    "time": "mean",
    "steiner_points": "mean",
    "energy": "mean",
    "p_bar": "mean",
    "obtuse_count": "mean"
}).reset_index()

##########################
# 2) Randomization Summary
##########################
randomization_summary = data.groupby("randomization").agg({
    "time": "mean",
    "steiner_points": "mean",
    "energy": "mean",
    "p_bar": "mean"
}).reset_index()

##########################
# 3) Overall Method Performance
##########################
methods = data["method"].unique()
method_performance = data.groupby("method").agg({
    "time": "mean",
    "steiner_points": "mean",
    "energy": "mean",
    "p_bar": "mean"
}).reset_index()

method_performance_energy = data.groupby("method").agg({
    "time": "mean",
    "steiner_points": "mean",
    "energy": "mean",
    "p_bar": "mean"
}).reset_index()

##########################
# 4) Heatmap (Method x Category)
##########################
heatmap_data = data.groupby(["method", "category"]).agg({
    "time": "mean",
    "steiner_points": "mean",
    "energy": "mean",
    "p_bar": "mean"
}).unstack().fillna(0)

##########################
# 5) Trend of p_bar by Category
##########################
category_p_bar_trend = data.groupby("category")["p_bar"].mean().reset_index()

##########################
# 6) Randomization Impact
##########################
randomization_impact = data.groupby("randomization").agg({
    "time": "mean",
    "steiner_points": "mean",
    "energy": "mean",
    "p_bar": "mean"
}).reset_index()

##########################
# 7) Create Output Dirs
##########################
plots_dir = "Project3/scripts/output/plots"
summaries_dir = "Project3/scripts/output/summaries"
os.makedirs(plots_dir, exist_ok=True)
os.makedirs(summaries_dir, exist_ok=True)

##########################
# PLOTTING
##########################

# 1. Average Time by Category
plt.figure(figsize=(8, 5))
sns.barplot(data=category_summary, x="category", y="time", hue="category", palette="coolwarm", legend=False)
plt.title("Average Time by Category")
plt.ylabel("Time (s)")
plt.xlabel("Category")
plt.savefig(os.path.join(plots_dir, "average_time_by_category.png"))
plt.show()


# 2. Performance of Randomized vs Non-Randomized Methods
plt.figure(figsize=(8, 5))
sns.barplot(data=randomization_summary, x="randomization", y="energy", hue="randomization", palette="coolwarm", legend=False)
plt.title("Average Energy: Randomized vs Non-Randomized Methods")
plt.ylabel("Energy")
plt.xlabel("Randomization")
plt.savefig(os.path.join(plots_dir, "average_energy_randomized_vs_non_randomized.png"))
plt.show()

# 3. Method Performance
plt.figure(figsize=(10, 6))
sns.barplot(data=method_performance, x="method", y="p_bar", hue="method", palette="coolwarm", legend=False)
plt.title("Performance (p̅) of Methods")
plt.ylabel("Average Convergence Rate (p̅)")
plt.xlabel("Method")
plt.savefig(os.path.join(plots_dir, "method_performance.png"))
plt.show()

# 4. Method Energy Performance
plt.figure(figsize=(10, 6))
sns.barplot(data=method_performance_energy, x="method", y="energy", hue="method", palette="coolwarm", legend=False)
plt.title("Energy Performance of Methods")
plt.ylabel("Average Energy")
plt.xlabel("Method")
plt.savefig(os.path.join(plots_dir, "method_energy_performance.png"))
plt.show()

# 5. Heatmap: Aggregated Performance Metrics by Method and Category
plt.figure(figsize=(12, 8))
sns.heatmap(heatmap_data["energy"], annot=True, fmt=".2f", cmap="coolwarm")
plt.title("Heatmap of Average Energy by Method and Category")
plt.xlabel("Category")
plt.ylabel("Method")
plt.savefig(os.path.join(plots_dir, "heatmap_energy_by_method_and_category.png"))
plt.show()

# 6. Trend of Average Convergence Rate Across Categories
plt.figure(figsize=(10, 6))
sns.lineplot(data=category_p_bar_trend, x="category", y="p_bar", marker="o", palette="coolwarm")
plt.title("Trend of Average Convergence Rate Across Categories")
plt.xlabel("Category")
plt.ylabel("Average Convergence Rate (p̅)")
plt.savefig(os.path.join(plots_dir, "trend_convergence_rate_by_category.png"))
plt.show()

# 7. Energy vs. Average Convergence Rate (p̅)
plt.figure(figsize=(8, 6))
sns.scatterplot(data=data, x="energy", y="p_bar", hue="method", style="randomization", palette="coolwarm")
plt.title("Energy vs. Average Convergence Rate (p̅)")
plt.xlabel("Energy")
plt.ylabel("Average Convergence Rate (p̅)")
plt.legend(title="Method & Randomization", loc="best")
plt.savefig(os.path.join(plots_dir, "energy_vs_convergence_rate.png"))
plt.show()

# 8. Effectiveness of Randomization: Energy vs. Average Convergence Rate (p̅)
plt.figure(figsize=(8, 6))
sns.scatterplot(data=data, x="energy", y="p_bar", hue="randomization", style="method", palette="coolwarm")
plt.title("Effectiveness of Randomization: Energy vs. Average Convergence Rate (p̅)")
plt.xlabel("Energy")
plt.ylabel("Average Convergence Rate (p̅)")
plt.legend(title="Randomization & Method", loc="best")
plt.savefig(os.path.join(plots_dir, "randomization_effectiveness_energy_vs_convergence.png"))
plt.show()

# 9. Impact of Randomization on Key Metrics
plt.figure(figsize=(8, 6))
sns.barplot(data=randomization_impact.melt(id_vars=["randomization"], var_name="Metric", value_name="Value"), 
            x="Metric", y="Value", hue="randomization", palette="coolwarm")
plt.title("Impact of Randomization on Key Metrics")
plt.ylabel("Average Value")
plt.xlabel("Performance Metric")
plt.xticks(rotation=45)
plt.savefig(os.path.join(plots_dir, "randomization_impact_metrics.png"))
plt.show()

# 10. Energy Distribution Across Randomized and Non-Randomized Methods
plt.figure(figsize=(10, 6))
sns.histplot(data=data, x="energy", hue="randomization", kde=True, palette="coolwarm", bins=20)
plt.title("Energy Distribution Across Randomized and Non-Randomized Methods")
plt.xlabel("Energy")
plt.ylabel("Frequency")
plt.savefig(os.path.join(plots_dir, "energy_distribution_randomization.png"))
plt.show()

# 11. Pairplot: Key Metrics by Method
pairplot_data = data[["method", "time", "energy", "p_bar", "steiner_points"]]
sns.pairplot(pairplot_data, hue="method", palette="coolwarm", diag_kind="kde")
plt.savefig(os.path.join(plots_dir, "pairplot_key_metrics_by_method.png"))
plt.show()

# 12. Average Steiner Points by Method and Randomization
steiner_points_summary = data.groupby(["method", "randomization"]).agg({"steiner_points": "mean"}).reset_index()

plt.figure(figsize=(10, 6))
sns.barplot(data=steiner_points_summary, x="method", y="steiner_points", hue="randomization", palette="coolwarm")
plt.title("Average Steiner Points by Method and Randomization")
plt.xlabel("Method")
plt.ylabel("Average Steiner Points")
plt.savefig(os.path.join(plots_dir, "average_steiner_points_by_method_randomization.png"))
plt.show()

# 13. Time Distribution by Category and Method
plt.figure(figsize=(12, 6))
sns.stripplot(data=data, x="category", y="time", hue="method", dodge=True, palette="coolwarm", alpha=0.7)
plt.title("Time Distribution by Category and Method")
plt.xlabel("Category")
plt.ylabel("Time (s)")
plt.legend(title="Method", loc="upper right")
plt.savefig(os.path.join(plots_dir, "time_distribution_by_category_method.png"))
plt.show()

# 14. Time Distribution by Method
plt.figure(figsize=(10, 6))
sns.histplot(data=data, x="time", hue="method", kde=True, bins=20, palette="coolwarm", element="step")
plt.title("Time Distribution by Method")
plt.xlabel("Time (s)")
plt.ylabel("Frequency")
plt.savefig(os.path.join(plots_dir, "time_distribution_by_method.png"))
plt.show()

# 15. Energy and Steiner Points Across Categories
category_trends = data.groupby("category").agg({"energy": "mean", "steiner_points": "mean"}).reset_index()
plt.figure(figsize=(12, 6))
sns.lineplot(data=category_trends, x="category", y="energy", marker="o", label="Energy", color="blue")
sns.lineplot(data=category_trends, x="category", y="steiner_points", marker="o", label="Steiner Points", color="red")
plt.title("Energy and Steiner Points Across Categories")
plt.xlabel("Category")
plt.ylabel("Average Value")
plt.legend(title="Metric")
plt.savefig(os.path.join(plots_dir, "energy_steiner_points_by_category.png"))
plt.show()

########################################
#   ADDITIONAL SUMMARIES & CSV EXPORTS
########################################

########################################################################
# A) BEST METHOD PER CATEGORY (BASIC) – as before, but only by method
########################################################################
best_method_per_category_list = []
for cat in data["category"].unique():
    subdf = data[data["category"] == cat]
    if len(subdf) > 0:
        # Primary: highest p_bar
        # Secondary: if tie, fewer steiner_points
        # Tertiary: if still tie, lower energy
        # We can define a custom sort key:
        subdf_sorted = subdf.sort_values(
            by=["p_bar", "steiner_points", "energy"],
            ascending=[False, True, True]  # Highest p_bar, then lowest steiner_points, then lowest energy
        )
        best_row = subdf_sorted.iloc[0]
        best_method_per_category_list.append({
            "category": cat,
            "best_method": best_row["method"],
            "best_steiner_method": best_row["steiner_methods"],
            "p_bar": best_row["p_bar"],
            "steiner_points": best_row["steiner_points"],
            "energy": best_row["energy"],
            "time": best_row["time"]
        })

best_method_per_category_df = pd.DataFrame(best_method_per_category_list)
best_method_per_category_df.to_csv("Project3/scripts/output/summaries/best_method_per_category.csv", index=False)

########################################################################
# B) BEST COMBINATION PER CATEGORY
#    Incorporating BOTH "method" + "steiner_methods" as a combined approach
########################################################################
# We can treat the combination of (method, steiner_methods) as a single "method+steiner" approach
data["method_steiner_combo"] = data["method"] + " + " + data["steiner_methods"].astype(str)

best_combo_per_category_list = []
for cat in data["category"].unique():
    subdf = data[data["category"] == cat]
    if len(subdf) > 0:
        # Sort using the same performance criteria:
        subdf_sorted = subdf.sort_values(
            by=["p_bar", "steiner_points", "energy"],
            ascending=[False, True, True]
        )
        best_row = subdf_sorted.iloc[0]
        best_combo_per_category_list.append({
            "category": cat,
            "best_method_steiner_combo": best_row["method_steiner_combo"],
            "p_bar": best_row["p_bar"],
            "steiner_points": best_row["steiner_points"],
            "energy": best_row["energy"],
            "time": best_row["time"],
        })

best_combo_per_category_df = pd.DataFrame(best_combo_per_category_list)
best_combo_per_category_df.to_csv("Project3/scripts/output/summaries/best_combo_per_category.csv", index=False)

########################################################################
# C) RANDOMIZATION IMPACT SUMMARY
########################################################################
r_impact_dict = {}
avg_energy_rand = data[data["randomization"] == True]["energy"].mean()
avg_energy_nonrand = data[data["randomization"] == False]["energy"].mean()
avg_pbar_rand = data[data["randomization"] == True]["p_bar"].mean()
avg_pbar_nonrand = data[data["randomization"] == False]["p_bar"].mean()

r_impact_dict["Energy_Improvement"] = round(avg_energy_nonrand - avg_energy_rand, 2)
r_impact_dict["p_bar_Change"] = round(avg_pbar_nonrand - avg_pbar_rand, 2)

randomization_impact_summary_df = pd.DataFrame([r_impact_dict])
randomization_impact_summary_df.to_csv("Project3/scripts/output/summaries/randomization_impact_summary.csv", index=False)

########################################################################
# D) METHOD + STEINER_METHOD SUMMARY (Pivoted)
#    Summarize average p_bar, steiner_points, energy for each combination
########################################################################
method_steiner_summary = data.groupby(["method", "steiner_methods"]).agg({
    "p_bar": "mean",
    "steiner_points": "mean",
    "energy": "mean",
    "time": "mean"
}).reset_index()

method_steiner_summary.to_csv("Project3/scripts/output/summaries/method_steiner_summary.csv", index=False)

########################################################################
# E) METHOD x CATEGORY SUMMARY
########################################################################
method_category_summary = data.groupby(["method", "category"]).agg({
    "p_bar": "mean",
    "energy": "mean",
    "steiner_points": "mean",
    "time": "mean"
}).reset_index()

method_category_summary.to_csv("Project3/scripts/output/summaries/method_category_summary.csv", index=False)

########################################################################
# F) Summaries Dictionary (optional)
########################################################################
summary = {
    "category_summary": category_summary,
    "randomization_summary": randomization_summary,
    "method_performance": method_performance
}

########################################################################
# G) Generate Conclusions
########################################################################
def generate_conclusions(data, randomization_summary):
    """
    Generate conclusions based on method performance and randomization impact.
    
    Since we have multiple methods and steiner methods, 
    we focus on 'method' for high-level insight.
    """
    # Identify the best-performing method based on average convergence rate
    best_method_row = data.loc[data["p_bar"].idxmax()]
    best_method = best_method_row["method"]
    best_method_p_bar = best_method_row["p_bar"]

    # Determine randomization's impact on energy
    randomized_energy = randomization_summary.loc[randomization_summary["randomization"] == True, "energy"].mean()
    non_randomized_energy = randomization_summary.loc[randomization_summary["randomization"] == False, "energy"].mean()
    if randomized_energy < non_randomized_energy:
        randomization_effect = "improves"
    else:
        randomization_effect = "does not improve"

    # Determine the method with the lowest average energy
    lowest_energy_row = data.loc[data["energy"].idxmin()]
    lowest_energy_method = lowest_energy_row["method"]
    lowest_energy = lowest_energy_row["energy"]

    # Find the method with the highest time efficiency (lowest time)
    fastest_method_row = data.loc[data["time"].idxmin()]
    fastest_method = fastest_method_row["method"]
    fastest_time = fastest_method_row["time"]

    # Summary of convergence rate
    average_convergence = data["p_bar"].mean()

    # Construct the conclusion string
    conclusion = (
        f"1. The best-performing method overall is '{best_method}' "
        f"with an average convergence rate (p̅) of {best_method_p_bar:.2f}.\n"
        f"2. Randomization generally {randomization_effect} performance, "
        f"with average energy for randomized methods being {randomized_energy:.2f} "
        f"compared to {non_randomized_energy:.2f} for non-randomized methods.\n"
        f"3. The method with the lowest average energy is '{lowest_energy_method}' "
        f"with an energy value of {lowest_energy:.2f}.\n"
        f"4. The fastest method is '{fastest_method}' with an average execution time "
        f"of {fastest_time:.4f} seconds.\n"
        f"5. The overall average convergence rate across all methods is {average_convergence:.2f}.\n"
        f"6. Randomization's impact is particularly significant in certain categories "
        f"or methods, as indicated in the visualizations.\n"
    )
    
    return conclusion

print("Summary and Conclusions")
print(generate_conclusions(method_performance, randomization_summary))

# Export the initial main summaries to CSV
category_summary.to_csv("Project3/scripts/output/summaries/category_summary.csv", index=False)
randomization_summary.to_csv("Project3/scripts/output/summaries/randomization_summary.csv", index=False)
method_performance.to_csv("Project3/scripts/output/summaries/method_performance.csv", index=False)