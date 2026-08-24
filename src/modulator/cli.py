from __future__ import annotations

import argparse
from pathlib import Path

from modulator.configuration import apply_overrides, load_config
from modulator.pipeline import ModulatorPipeline, STAGE_ORDER
from modulator.runtime import find_project_root


def _parse_stage_list(raw: str | None) -> list[str] | None:
    if not raw or raw.strip().lower() == "all":
        return None
    stages = [item.strip() for item in raw.split(",") if item.strip()]
    bad = [stage for stage in stages if stage not in STAGE_ORDER]
    if bad:
        raise ValueError(f"Unknown stage(s): {', '.join(bad)}. Valid stages: {', '.join(STAGE_ORDER)}")
    return stages


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Package-first runner for the modulator pipeline.")
    sub = parser.add_subparsers(dest="command", required=True)

    run = sub.add_parser("run", help="Run the pipeline from a config file.")
    run.add_argument("--config", default="config/config.yaml", help="Path to YAML config relative to the project root or absolute.")
    run.add_argument("--workdir", default=".", help="Project directory containing workflow/, config/, results/, and resources/.")
    run.add_argument("--set", dest="overrides", nargs="*", default=[], help="Simple overrides like key=value or nested.key=value.")
    run.add_argument("--jobs", type=int, default=1, help="Number of independent sample-level jobs to run in parallel.")
    run.add_argument("--stages", default="all", help=f"Comma-separated subset of stages to run. Valid: {', '.join(STAGE_ORDER)}")
    run.add_argument("--resume", action="store_true", help="Skip stages whose outputs already exist in the results folder (checkpoint resume).")

    validate = sub.add_parser("validate-config", help="Load the config, apply overrides, and print the resolved project root.")
    validate.add_argument("--config", default="config/config.yaml")
    validate.add_argument("--workdir", default=".")
    validate.add_argument("--set", dest="overrides", nargs="*", default=[])

    demo = sub.add_parser("demo", help="Run a fast bundled demo dataset with explicit reference inputs.")
    demo.add_argument("--config", default="config/config.yaml", help="Base YAML config to start from.")
    demo.add_argument("--workdir", default=".", help="Project directory containing workflow/, config/, results/, and resources/.")
    demo.add_argument("--reference-fa", required=True, help="Reference FASTA to use for the demo run.")
    demo.add_argument("--reference-gtf", required=True, help="Reference GTF to use for the demo run.")
    demo.add_argument(
        "--dataset",
        # Only datasets actually bundled under resources/test_bams/ -- RPL13_reads was advertised but
        # never shipped, so `demo --dataset RPL13_reads` failed at bams_dir resolution.
        choices=["MXD1_reads", "ALCAM_NHSL1_SERAC1_MXD1_RIOK3_reads"],
        default="MXD1_reads",
        help="Bundled test dataset to use.",
    )
    demo.add_argument(
        "--mode",
        choices=["quick", "full"],
        default="quick",
        help="quick keeps the demo lightweight; full runs the standard stage set.",
    )
    demo.add_argument("--prefix", default="", help="Optional output prefix override.")
    demo.add_argument("--jobs", type=int, default=2, help="Number of independent sample-level jobs to run in parallel.")
    demo.add_argument("--set", dest="overrides", nargs="*", default=[], help="Additional key=value overrides.")
    return parser


def cmd_run(args: argparse.Namespace) -> None:
    project_root = find_project_root(args.workdir)
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = project_root / config_path
    config = apply_overrides(load_config(config_path), args.overrides)
    stages = _parse_stage_list(args.stages)
    pipeline = ModulatorPipeline(config, workdir=project_root, jobs=args.jobs, verbose=True, resume=args.resume)
    pipeline.run(stages=stages)


def cmd_validate(args: argparse.Namespace) -> None:
    project_root = find_project_root(args.workdir)
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = project_root / config_path
    config = apply_overrides(load_config(config_path), args.overrides)
    # validate-config must not mutate the filesystem: stage_inputs=False resolves + validates the
    # samplesheet/contrasts without staging BAM symlinks or writing the metadata TSV.
    pipeline = ModulatorPipeline(config, workdir=project_root, jobs=1, verbose=False, stage_inputs=False)
    print(f"project_root={pipeline.root}")
    print(f"prefix={pipeline.prefix}")
    print(f"samples={len(pipeline.samples)}")
    print(f"reference_fa={pipeline.reference_fa}")
    print(f"reference_gtf={pipeline.reference_gtf}")


def cmd_demo(args: argparse.Namespace) -> None:
    project_root = find_project_root(args.workdir)
    config_path = Path(args.config)
    if not config_path.is_absolute():
        config_path = project_root / config_path

    default_prefix = f"demo_{args.dataset.lower()}_{args.mode}"
    overrides = [
        f"reference_fa={args.reference_fa}",
        f"reference_gtf={args.reference_gtf}",
        f"bams_dir=resources/test_bams/{args.dataset}",
        f"prefix={args.prefix or default_prefix}",
        "assembler.write_zt_tagged_sample_bams=true",
        "genotype.enable=false",
    ]
    if args.mode == "quick":
        overrides.extend([
            "toggles.enable_zn_pileup=true",
            "toggles.enable_zn_aggregate=true",
            "toggles.enable_zt_pileup=false",
            "toggles.enable_zt_aggregate=false",
            "toggles.enable_test_diffs=false",
            "report.enable=true",
        ])
    config = apply_overrides(load_config(config_path), [*overrides, *args.overrides])
    pipeline = ModulatorPipeline(config, workdir=project_root, jobs=args.jobs, verbose=True)
    pipeline.run()


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()
    if args.command == "run":
        cmd_run(args)
    elif args.command == "validate-config":
        cmd_validate(args)
    elif args.command == "demo":
        cmd_demo(args)
    else:
        parser.error(f"Unhandled command: {args.command}")


if __name__ == "__main__":
    main()
