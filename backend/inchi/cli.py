import argparse
import json

from backend.inchi.compare import compare_text_files, compare_pair, compare_mgf_files
from backend.inchi.config_loader import load_config, build_config_from_levels, apply_inchitrust

def main():
    parser = argparse.ArgumentParser(
        prog="inchi",
        description="Metabolite identity comparison tool"
    )

    subparsers = parser.add_subparsers(dest="command")

    compare_parser = subparsers.add_parser(
        "compare",
        help="Compare two txt files containing InChIs"
    )

    compare_parser.add_argument("file1")
    compare_parser.add_argument("file2")

    compare_parser.add_argument(
        "--config",
        help="Path to config.json",
        default=None
    )

    compare_parser.add_argument(
        "--output_file",
        help="Output results JSON",
        default=None
    )

    add_inchitrust_arg(compare_parser)

    compare_parser.add_argument(
        "--only-differences",
        action="store_true",
        help="Show only comparisons with differences"
    )

    pair_parser = subparsers.add_parser(
        "compare-pair",
        help="Compare two InChIs"
    )

    pair_parser.add_argument("inchi1")
    pair_parser.add_argument("inchi2")

    pair_parser.add_argument(
        "--config",
        help="Path to config.json",
        default=None
    )

    add_inchitrust_arg(pair_parser)

    pair_levels_parser = subparsers.add_parser(
        "compare-pair-levels",
        help="Compare two InChIs with selected identity levels"
    )

    pair_levels_parser.add_argument("inchi1")
    pair_levels_parser.add_argument("inchi2")

    pair_levels_parser.add_argument(
        "--levels",
        nargs="+",
        required=True,
        help="List of levels (e.g. complete_identity isotope charge)"
    )

    pair_levels_parser.add_argument(
        "--config",
        help="Base config.json",
        default=None
    )

    add_inchitrust_arg(pair_levels_parser)

    mgf_parser = subparsers.add_parser(
        "compare-mgf",
        help="Compare two MGF metabolomics files"
    )

    mgf_parser.add_argument("file1")
    mgf_parser.add_argument("file2")

    mgf_parser.add_argument(
        "--config",
        help="Path to config.json",
        default=None
    )

    add_inchitrust_arg(mgf_parser)

    args = parser.parse_args()

    if args.command == "compare":
        config = load_config(args.config)

        result = compare_text_files(
            args.file1,
            args.file2,
            config,
            only_differences=args.only_differences
        )

        if args.output_file:
            with open(args.output_file, "w") as f:
                json.dump(result, f, indent=2)
        else:
            print(json.dumps(result, indent=2))

    elif args.command == "compare-pair":
        config = load_config(args.config)
        config = apply_inchitrust(config, args.inchitrust)

        result = compare_pair(
            args.inchi1,
            args.inchi2,
            config
        )

        print(json.dumps(result, indent=2))

    elif args.command == "compare-pair-levels":
        base_config = load_config(args.config)
        base_config = apply_inchitrust(config, args.inchitrust)

        config = build_config_from_levels(
            args.levels,
            base_config
        )

        result = compare_pair(
            args.inchi1,
            args.inchi2,
            config
        )

        print(json.dumps(result, indent=2))

    elif args.command == "compare-mgf":
        config = load_config(args.config)
        config = apply_inchitrust(config, args.inchitrust)

        result = compare_mgf_files(
            args.file1,
            args.file2,
            config
        )

        print(json.dumps(result, indent=2))

    else:
        parser.print_help()

def add_inchitrust_arg(parser):
    parser.add_argument(
        "--inchitrust",
        help="Path to InChI Trust executable",
        default=None
    )

if __name__ == "__main__":
    main()