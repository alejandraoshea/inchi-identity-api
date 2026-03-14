import argparse
from inchi.compare import compare_files

def main():
    parser = argparse.ArgumentParser(
        prog="inchi",
        description="Inchi identity comparison tool"
    )

    subparsers = parser.add_subparsers(dest="command")

    # compare command
    compare_parser = subparsers.add_parser(
        "compare",
        help="Compare two files containing InChIs"
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
        default="output/results.json"
    )

    args = parser.parse_args()

    if args.command == "compare":
        compare_files(
            args.file1,
            args.file2,
            args.config,
            args.output_file
        )

if __name__ == "__main__":
    main()