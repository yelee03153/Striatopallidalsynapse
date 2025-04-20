import seaborn as sns
import matplotlib.pyplot as plt
from scipy.stats import mannwhitneyu, wilcoxon

# Parameters to analyze
parameters = ['calcium_channel_open_prob', 'total_Ca']

def plot_parameter_comparisons(df_all, region, exp_type):
    fig, axes = plt.subplots(len(parameters), 1, figsize=(10, len(parameters)*4))
    region_data = df_all[df_all['Region'] == region]

    for i, param in enumerate(parameters):
        sns.violinplot(data=region_data, x='Condition', y=param, ax=axes[i], inner=None, palette="Set2")
        sns.scatterplot(data=region_data, x='Condition', y=param, ax=axes[i], color="black", s=50)
        axes[i].set_title(f'{param} - {region}')

        control_data = region_data[region_data['Condition'] == 'Baseline'][param]
        quinpirole_data = region_data[region_data['Condition'] == 'Quinpirole'][param]

        stat, p_value = mannwhitneyu(control_data, quinpirole_data)
        significance = ('***' if p_value < 0.001 else
                        '**' if p_value < 0.01 else
                        '*' if p_value < 0.05 else 'ns')

        axes[i].text(0.5, 1.05, f'p = {p_value:.4f} {significance}', ha='center', va='bottom',
                     transform=axes[i].transAxes, fontsize=12)

    plt.tight_layout()
    plt.savefig(f"results/Plots/Stat_{region}_{exp_type}.pdf")
    plt.show()


def plot_paired_comparisons(df_all, region,exp_type):
    fig, axes = plt.subplots(len(parameters), 1, figsize=(10, len(parameters)*4))
    region_data = df_all[df_all['Region'] == region]

    for i, param in enumerate(parameters):
        control_data = region_data[region_data['Condition'] == 'Baseline'][param]
        quinpirole_data = region_data[region_data['Condition'] == 'Quinpirole'][param]

        sns.violinplot(data=region_data[region_data['Condition'] == 'Baseline'], x='Condition', y=param,
                       ax=axes[i], inner=None, color='gray', alpha=0.3)
        sns.violinplot(data=region_data[region_data['Condition'] == 'Quinpirole'], x='Condition', y=param,
                       ax=axes[i], inner=None, color='red', alpha=0.3)

        for j in range(len(control_data)):
            axes[i].plot([0, 1], [control_data.iloc[j], quinpirole_data.iloc[j]], color='red', alpha=0.5)

        try:
            stat, p_value = wilcoxon(control_data, quinpirole_data, zero_method='zsplit')
        except ValueError as e:
            if "x - y is zero for all elements" in str(e):
                p_value = 1.0
            else:
                raise e

        significance = ('***' if p_value < 0.001 else
                        '**' if p_value < 0.01 else
                        '*' if p_value < 0.05 else 'ns')

        axes[i].text(0.5, 1.05, f'p = {p_value:.4f} {significance}', ha='center', va='bottom',
                     transform=axes[i].transAxes, fontsize=12)

        axes[i].set_ylabel("Measurement Value")

    plt.tight_layout()
    plt.savefig(f"results/Plots/Stat_{region}_pairwis_{exp_type}.pdf")
    plt.show()
