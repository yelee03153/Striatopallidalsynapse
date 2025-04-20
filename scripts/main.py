import argparse
from analysis import run_fitting_pipeline
from vis import plot_parameter_comparisons

def main():
    parser = argparse.ArgumentParser(description="Run synaptic fitting and analysis pipeline.")
    parser.add_argument("--region", type=str, default="VM", help="Brain region to analyze (e.g., VM)")
    parser.add_argument("--n", type=int, default=10, help="Number of valid fits to collect")
    parser.add_argument("--exp_type", type=str, default="control", help="Type of experiment (control or 6OHDA)")

    args = parser.parse_args()
    region = args.region
    n = args.n
    exp_type=args.exp_type
    # Run analysis
    df_all = run_fitting_pipeline(region, n)

    # Save parameter DataFrame
    df_all.to_csv(f"results/Output/DF_All_{region}_{exp_type}.csv", index=False)

    # Plot results
    plot_parameter_comparisons(df_all, region, exp_type)

if __name__ == "__main__":
    main()
