#!/usr/bin/env python3
"""Create a Steinbock-compatible panel.csv from raw IMC ROI text files.

The script scans ROI `.txt` files in a raw data directory, extracts channel and
marker names from the first header, validates that all ROI files share the same
channel order, and writes a dynamic Steinbock panel.

Default DeepCell/Mesmer assignments for this workflow are:

- nuclear markers/channels: HistoneH3, 191Ir, 193Ir
- membrane markers: CD98, CD3, CD138, CD45

These defaults can be changed with command-line arguments when applying the
workflow to a new IMC dataset.
"""

from __future__ import annotations

import argparse
import csv
import re
import shutil
from pathlib import Path


METADATA_COLUMNS = {"Start_push", "End_push", "Pushes_duration", "X", "Y", "Z"}
CHANNEL_PATTERN = re.compile(r"^(?P<marker>.+)\((?P<channel>.+?)(?:Di|Ir)\)$")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Create data/panel.csv from raw IMC ROI text-file headers."
    )
    parser.add_argument(
        "--workflow-dir",
        type=Path,
        default=Path.cwd(),
        help="Workflow directory containing data/, results/, scripts/, and notebooks/.",
    )
    parser.add_argument(
        "--raw-dir",
        type=Path,
        default=Path("data/raw"),
        help="Directory containing raw .mcd files and ROI .txt files.",
    )
    parser.add_argument(
        "--output",
        type=Path,
        default=Path("data/panel.csv"),
        help="Output Steinbock panel CSV path. Steinbock expects data/panel.csv.",
    )
    parser.add_argument(
        "--numbered-copy",
        type=Path,
        default=Path("results/01_panel.csv"),
        help=(
            "Chronologically numbered research copy. "
            "Use 'none' to skip creating a numbered copy."
        ),
    )
    parser.add_argument(
        "--nuclear-marker",
        nargs="+",
        default=["HistoneH3", "191Ir", "193Ir"],
        help="Marker(s) or channel label(s) assigned as DeepCell/Mesmer channel 1.",
    )
    parser.add_argument(
        "--membrane-marker",
        nargs="+",
        default=["CD98", "CD3", "CD138", "CD45"],
        help="Marker(s) or channel label(s) assigned as DeepCell/Mesmer channel 2.",
    )
    return parser.parse_args()


def workflow_path(workflow_dir: Path, path: Path) -> Path:
    return path if path.is_absolute() else workflow_dir / path


def read_header(path: Path) -> list[str]:
    with path.open("r", encoding="utf-8", errors="replace") as handle:
        return handle.readline().rstrip("\n").split("\t")


def extract_marker_channels(header: list[str]) -> list[tuple[str, str]]:
    marker_channels: list[tuple[str, str]] = []

    for column in header:
        if column in METADATA_COLUMNS:
            continue

        match = CHANNEL_PATTERN.match(column)
        if match is None:
            continue

        marker = match.group("marker").strip()
        channel = match.group("channel").strip()
        marker_channels.append((channel, marker))

    if not marker_channels:
        raise ValueError("No IMC marker/channel columns were detected in the ROI header.")

    return marker_channels


def find_roi_text_files(raw_dir: Path) -> list[Path]:
    roi_files = sorted(raw_dir.glob("*.txt"))
    if not roi_files:
        raise FileNotFoundError(f"No ROI .txt files were found in {raw_dir}")
    return roi_files


def validate_consistent_headers(roi_files: list[Path]) -> list[tuple[str, str]]:
    reference = extract_marker_channels(read_header(roi_files[0]))
    mismatches: list[str] = []

    for roi_file in roi_files[1:]:
        current = extract_marker_channels(read_header(roi_file))
        if current != reference:
            mismatches.append(str(roi_file))

    if mismatches:
        mismatch_list = "\n".join(mismatches)
        raise ValueError(
            "Not all ROI files share the same marker/channel order. "
            f"Mismatching files:\n{mismatch_list}"
        )

    return reference


def deepcell_value(
    marker: str,
    channel: str,
    nuclear_markers: list[str],
    membrane_markers: list[str],
) -> str:
    if marker in nuclear_markers or channel in nuclear_markers:
        return "1"
    if marker in membrane_markers or channel in membrane_markers:
        return "2"
    return ""


def write_panel(
    marker_channels: list[tuple[str, str]],
    output: Path,
    nuclear_markers: list[str],
    membrane_markers: list[str],
) -> None:
    output.parent.mkdir(parents=True, exist_ok=True)

    with output.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=["channel", "name", "keep", "ilastik", "deepcell", "cellpose"],
        )
        writer.writeheader()

        for index, (channel, marker) in enumerate(marker_channels, start=1):
            writer.writerow(
                {
                    "channel": channel,
                    "name": marker,
                    "keep": "1",
                    "ilastik": str(index),
                    "deepcell": deepcell_value(
                        marker,
                        channel,
                        nuclear_markers,
                        membrane_markers,
                    ),
                    "cellpose": "",
                }
            )


def main() -> None:
    args = parse_args()
    workflow_dir = args.workflow_dir.resolve()
    raw_dir = workflow_path(workflow_dir, args.raw_dir).resolve()
    output = workflow_path(workflow_dir, args.output).resolve()
    numbered_copy = (
        None
        if str(args.numbered_copy).casefold() == "none"
        else workflow_path(workflow_dir, args.numbered_copy).resolve()
    )

    roi_files = find_roi_text_files(raw_dir)
    marker_channels = validate_consistent_headers(roi_files)
    write_panel(marker_channels, output, args.nuclear_marker, args.membrane_marker)

    if numbered_copy is not None:
        numbered_copy.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(output, numbered_copy)

    print(f"ROI text files checked: {len(roi_files)}")
    print(f"Panel rows written: {len(marker_channels)}")
    print(f"Steinbock panel path: {output}")
    if numbered_copy is not None:
        print(f"Numbered research copy: {numbered_copy}")
    print(f"DeepCell channel 1 markers/channels: {args.nuclear_marker}")
    print(f"DeepCell channel 2 markers/channels: {args.membrane_marker}")


if __name__ == "__main__":
    main()
