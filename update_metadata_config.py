import json
import argparse
import sys


def option_parser():
    parser = argparse.ArgumentParser(description="update metadata config")
    parser.add_argument('--config', "-m",
                        help='metadata json config file', dest="config")
    parser.add_argument('--condition_1', "-a",
                        help='condition 1, format e.g. name@sample', dest="condition_1", 
                        default="default", nargs='?')
    parser.add_argument('--condition_2', "-b",
                        help='condition 2, format e.g. name@sample', dest="condition_2", default="default",
                        nargs='?')
    parser.add_argument('--condition_3', "-c",
                        help='condition 3, format e.g. name@sample', dest="condition_3", default="default",
                        nargs='?')

    args = parser.parse_args()

    if len(sys.argv) < 2:
        parser.print_help()
        sys.exit("Missing parameters!")

    return args


def get_field(condition, header_name):
    tmp = condition.split("@")
    if len(tmp) == 2:
        object_type = tmp[1]
        field_name = tmp[0]
        field = [
            field_name
        ]

        if field_name == "timepoint":
            field = [
                "timepoints",
                "$TIMEPOINT",
                "timepoint_termId"
            ]
        return {
            "field": field,
            "object_type": object_type,
            "header_name": header_name
        }
    else:
        raise ValueError(
            "Condition needs to have the format field_name@object_type")


def main():
    args = option_parser()
    json_dict = json.load(open(args.config))
    columns = json_dict.get("columns")

    if args.condition_1 != "default":
        columns[3] = get_field(args.condition_1, columns[3].get("header_name"))
    if args.condition_2 != "default":
        columns[4] = get_field(args.condition_2, columns[4].get("header_name"))
    if args.condition_3 != "default":
        columns[5] = get_field(args.condition_3, columns[5].get("header_name"))
    json_dict["columns"] = columns
    json.dump(json_dict, open("meta_table_config_changed.json", "w"), indent=2)


if __name__ == "__main__":
    main()
