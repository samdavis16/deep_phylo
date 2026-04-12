import argparse
import os

from . import workflows

"""
NOTE: This interface calls a multi-iteration run from workflows.run_n_iters with only limited exposure of functionality 
(mainly uses default parameters / tools at each stage of the workflow, but can be applied to any input profile/s). A 
more comprehensive CLI is under development. Some aspects of the workflow are not yet optimised. Further options for 
some tools are still to be fully integrated (e.g. only MAFFT is currently available as an aligner).
"""


def main():

    parser = argparse.ArgumentParser()

    parser.add_argument(
        "--name",
        required=True,
        help="Name of run, used as prefix for output files.",
        metavar="NAME"
    )

    parser.add_argument(
        "--hmm",
        nargs=2,
        required=True,
        action="append",
        metavar=("HMM_FILE", "THRESHOLD"),
        help="Input profile files and thresholds. "
             "Specify as: --hmm <HMM_FILE> <THRESHOLD> for each profile."
    )

    parser.add_argument(
        "--iters",
        type=int,
        default=5
    )

    parser.add_argument(
        "--conv-prop",
        type=float,
        help="Defines convergence of exploration where "
             "new sequences as  of total < CONV-PROP. "
             "[0.0-1.0]. --iters is ignored if specified."
    )

    parser.add_argument(
        "--max-iters",
        type=int,
        default=10,
        help="Maximum number of iterations to run if convergence "
             "proportion is specified."
    )

    parser.add_argument(
        "--db-dir",
        help="Relative path to dir containing database for curation. "
             "Note: only Uniparc database currently supported. "
             "Downloads all Uniparc segments if not specified."
    )

    parser.add_argument(
        "--cpu",
        type=int,
        default=os.cpu_count(),
        help="# CPUs to employ. "
             "Maximum available on system by default."
    )

    args = parser.parse_args()

    # Unzip profiles / thresholds
    init_profiles = [prof_data[0] for prof_data in args.hmm]
    init_thresholds = [float(prof_data[1]) for prof_data in args.hmm]

    if args.conv_prop:
        args.iters = None

    cpu_per_search = min(6, args.cpu)
    max_search_processes = args.cpu // cpu_per_search

    workflows.run_n_iters(
        init_profiles,
        init_thresholds,
        args.name,
        n_full_iters=args.iters,
        convergence_prop=args.conv_prop,
        max_iters=args.max_iters,
        full_db_dir=args.db_dir,
        segment_idx_suffix=".faidx",
        max_search_processes=max_search_processes,
        cpu_per_search=cpu_per_search,
)

if __name__ == "__main__":
    main()