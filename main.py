import argparse

from config.paths import load_paths
from config.pipeline_config import load_pipeline_config
from src.workflow.runners import run_pipeline


def main():
    parser = argparse.ArgumentParser(description="PHDeep2.0 entrypoint")
    parser.add_argument(
        "--show-paths",
        action="store_true",
        help="Print resolved runtime paths and exit.",
    )
    parser.add_argument(
        "--validate-config",
        type=str,
        help="Validate a pipeline JSON config and print normalized summary.",
    )
    parser.add_argument(
        "--run-config",
        type=str,
        help="Run pipeline stages from a pipeline JSON config.",
    )
    args = parser.parse_args()

    if args.show_paths:
        paths = load_paths()
        print(f"data_root={paths.data_root}")
        print(f"output_root={paths.output_root}")
        print(f"reference_root={paths.reference_root}")
        print(f"logs_root={paths.logs_root}")
        return

    if args.validate_config:
        config = load_pipeline_config(args.validate_config)
        print(f"config_version={config.version}")
        print(f"backend={config.execution.backend}")
        print(f"stage_order={','.join(config.stage_order)}")
        return

    if args.run_config:
        config = load_pipeline_config(args.run_config)
        results = run_pipeline(config)
        print(f"backend={config.execution.backend}")
        for stage in config.stage_order:
            print(f"stage={stage}")
            for key, value in results.get(stage, {}).items():
                print(f"{stage}.{key}={value}")
        return

    print("PHDeep2.0 is installed. Use module CLIs for specific tasks.")


if __name__ == "__main__":
    main()
