import argparse
import json

from backend.inchi.compare import compare_text_files,compare_pair,compare_mgf_files
from backend.inchi.config_loader import load_config,build_config_from_levels,apply_inchitrust

def main():
    parser = argparse.ArgumentParser(
        prog="inchi",
        description="Metabolite identity comparison tool"
    )

    subparsers = parser.add_subparsers(dest="command")

    compare_parser = subparsers.add_parser("compare")

    compare_parser.add_argument("file1")
    compare_parser.add_argument("file2")

    compare_parser.add_argument(
        "--mode",
        choices=["pairwise", "cross"],
        default="pairwise"
    )

    compare_parser.add_argument("--config", default=None)

    compare_parser.add_argument(
        "--only-equal",
        action="store_true"
    )

    compare_parser.add_argument("--output_file", default=None)

    add_inchitrust_arg(compare_parser)

    pair_parser = subparsers.add_parser("compare-pair")

    pair_parser.add_argument("inchi1")
    pair_parser.add_argument("inchi2")
    pair_parser.add_argument("--config", default=None)

    add_inchitrust_arg(pair_parser)

    pair_levels_parser = subparsers.add_parser("compare-pair-levels")

    pair_levels_parser.add_argument("inchi1")
    pair_levels_parser.add_argument("inchi2")

    pair_levels_parser.add_argument(
        "--levels",
        nargs="+",
        required=True
    )

    pair_levels_parser.add_argument("--config", default=None)
    add_inchitrust_arg(pair_levels_parser)

    mgf_parser = subparsers.add_parser("compare-mgf")

    mgf_parser.add_argument("file1")
    mgf_parser.add_argument("file2")
    mgf_parser.add_argument("--config", default=None)

    mgf_parser.add_argument(
        "--only-equal",
        action="store_true"
    )
    mgf_parser.add_argument(
        "--level",
        default="COMPLETE_IDENTITY"
    )

    add_inchitrust_arg(mgf_parser)
    args = parser.parse_args()

    if args.command == "compare":
        config = load_config(args.config)
        config = apply_inchitrust(config, args.inchitrust)

        result = compare_text_files(
            args.file1,
            args.file2,
            config,
            mode=args.mode,
            only_equal=args.only_equal
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
        base_config = apply_inchitrust(base_config, args.inchitrust)

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
            config,
            level=args.level
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