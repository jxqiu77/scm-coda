#!/usr/bin/env python3
"""Convert emp_theo Excel output to LaTeX table in target paper style.

Usage:
  python3 data/xlsx_to_latex.py
  python3 data/xlsx_to_latex.py --dataset spike
"""

from __future__ import annotations

import argparse
import zipfile
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Dict, List, Tuple
from xml.etree import ElementTree as ET

NS = {"a": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


@dataclass
class DataRow:
    p: str
    n_total: str
    empirical: List[str]
    theoretical: List[str]
    dist: str
    c: str


@dataclass
class TableConfig:
    xlsx: str
    out: str
    caption: str
    label: str
    header_line1: str
    header_line2: str
    header_line3: str


TABLE_CONFIGS: Dict[str, TableConfig] = {
    "bidiagonal": TableConfig(
        xlsx="emp_theo_Bidiagonal.xlsx",
        out="emp_theo_Bidiagonal_table.tex",
        caption=r"\caption{Empirical size and power ($2000$ replications) for $\bSigma_2$; values in parentheses are based on known $\{w_{ij}\}$ distribution and true parameters are used.}",
        label=r"\label{tab:covtest-Sigma2}",
        header_line1=r"&   &     &     & \multicolumn{1}{c}{Size} & \multicolumn{1}{c}{} & \multicolumn{4}{c}{Power} \\",
        header_line2=r"\cmidrule(lr){5-5} \cmidrule(l){7-10}",
        header_line3=r"$c$ & $p$ & $n_1$ & $n_2$ & \multicolumn{1}{c}{{\tiny $\psi_2=0$}} & \multicolumn{1}{c}{} & \multicolumn{1}{c}{{\tiny $\psi_2=-0.2$}} & \multicolumn{1}{c}{{\tiny $\psi_2=-0.25$}} & \multicolumn{1}{c}{{\tiny $\psi_2=-0.3$}} & \multicolumn{1}{c}{{\tiny $\psi_2=-0.35$}} \\",
    ),
    "spike": TableConfig(
        xlsx="emp_theo_Spike.xlsx",
        out="emp_theo_Spike_table.tex",
        caption=r"\caption{Empirical size and power ($2000$ replications) for $\bSigma_{1}$; values in parentheses are based on known $\{w_{ij}\}$ distribution and true parameters are used.}",
        label=r"\label{tab:covtest-Sigma1}",
        header_line1=r"& \multicolumn{1}{c}{} & \multicolumn{1}{c}{} & \multicolumn{1}{c}{} & \multicolumn{1}{c}{Size} & \multicolumn{1}{c}{} & \multicolumn{4}{c}{Power} \\",
        header_line2=r"\cmidrule(lr){5-5} \cmidrule(l){7-10}",
        header_line3=r"$c$ & \multicolumn{1}{c}{$p$} & \multicolumn{1}{c}{$n_{1}$} & \multicolumn{1}{c}{$n_{2}$} & \multicolumn{1}{c}{$\psi_{1}=1$} & \multicolumn{1}{c}{} & \multicolumn{1}{c}{$\psi_{1}=3$} & \multicolumn{1}{c}{$\psi_{1}=3.5$} & \multicolumn{1}{c}{$\psi_{1}=4$} & \multicolumn{1}{c}{$\psi_{1}=4.5$} \\",
    ),
}


def _col_index(cell_ref: str) -> int:
    letters = "".join(ch for ch in cell_ref if ch.isalpha())
    idx = 0
    for ch in letters:
        idx = idx * 26 + (ord(ch.upper()) - 64)
    return idx


def _read_sheet_rows(xlsx_path: Path) -> List[List[str]]:
    with zipfile.ZipFile(xlsx_path) as zf:
        shared_strings: List[str] = []
        if "xl/sharedStrings.xml" in zf.namelist():
            sst = ET.fromstring(zf.read("xl/sharedStrings.xml"))
            for si in sst.findall("a:si", NS):
                text = "".join((t.text or "") for t in si.findall(".//a:t", NS))
                shared_strings.append(text)

        sheet_xml = ET.fromstring(zf.read("xl/worksheets/sheet1.xml"))
        rows: List[List[str]] = []
        for row in sheet_xml.findall(".//a:sheetData/a:row", NS):
            values: Dict[int, str] = {}
            for cell in row.findall("a:c", NS):
                cell_type = cell.attrib.get("t")
                value_node = cell.find("a:v", NS)
                value = ""
                if value_node is not None and value_node.text is not None:
                    value = value_node.text
                    if cell_type == "s":
                        value = shared_strings[int(value)]
                values[_col_index(cell.attrib["r"])] = value
            max_col = max(values) if values else 0
            rows.append([values.get(i, "") for i in range(1, max_col + 1)])
    return rows


def _to_int_str_if_possible(v: str) -> str:
    try:
        f = float(v)
    except ValueError:
        return v
    if f.is_integer():
        return str(int(f))
    return v


def _format_c_value(v: str) -> str:
    try:
        f = float(v)
    except ValueError:
        return v
    return str(int(f)) if f.is_integer() else str(f)


def parse_rows_to_records(rows: List[List[str]]) -> Tuple[List[str], List[DataRow]]:
    if not rows:
        raise ValueError("Excel is empty.")
    header = rows[0]
    if len(header) < 9:
        raise ValueError("Unexpected Excel format. Expected at least 9 columns.")

    records: List[DataRow] = []
    i = 1
    while i < len(rows):
        row = rows[i]
        if len(row) < 9:
            i += 1
            continue
        if row[0].strip():
            empirical = row[2:7]
            dist = row[7].strip()
            c_val = _format_c_value(row[8].strip())
            p_val = _to_int_str_if_possible(row[0].strip())
            n_val = _to_int_str_if_possible(row[1].strip())
            theoretical = [""] * 5

            if i + 1 < len(rows):
                next_row = rows[i + 1]
                if len(next_row) >= 7 and not next_row[0].strip():
                    theoretical = next_row[2:7]
                    i += 1

            records.append(
                DataRow(
                    p=p_val,
                    n_total=n_val,
                    empirical=empirical,
                    theoretical=theoretical,
                    dist=dist,
                    c=c_val,
                )
            )
        i += 1
    return header, records


def _sort_key_num(s: str) -> float:
    try:
        return float(s)
    except ValueError:
        return float("inf")


def build_table_body(records: List[DataRow]) -> List[str]:
    # Keep readable model names in output.
    model_label = {
        "Exp": "Exp(5)",
        "ChiSq": r"$\chi^2(1)$",
    }

    grouped: Dict[str, Dict[str, List[DataRow]]] = defaultdict(lambda: defaultdict(list))
    for rec in records:
        grouped[rec.dist][rec.c].append(rec)

    dist_order = sorted(grouped.keys(), key=lambda d: (0 if d == "Exp" else 1, d))
    lines: List[str] = []

    for dist_idx, dist in enumerate(dist_order):
        pretty_dist = model_label.get(dist, dist)
        lines.append(f"        &   &     &     & \\multicolumn{{6}}{{l}}{{{pretty_dist} model}} \\\\")

        c_values = sorted(grouped[dist].keys(), key=_sort_key_num)
        for c_idx, c in enumerate(c_values):
            rows_c = sorted(grouped[dist][c], key=lambda r: _sort_key_num(r.p))
            lines.append(f"        \\multirow{{{2 * len(rows_c)}}}{{*}}{{{c}}}")

            for rec in rows_c:
                n1 = n2 = ""
                try:
                    n_total = int(float(rec.n_total))
                    n1 = str(n_total // 2)
                    n2 = str(n_total // 2)
                except ValueError:
                    pass

                lines.append(
                    "            "
                    f"& {rec.p} & {n1} & {n2} & {rec.empirical[0]} &       & "
                    f"{rec.empirical[1]} & {rec.empirical[2]} & {rec.empirical[3]} & {rec.empirical[4]} \\\\"
                )
                lines.append(
                    "            "
                    f"&     &     &     & {rec.theoretical[0]} &       & "
                    f"{rec.theoretical[1]} & {rec.theoretical[2]} & {rec.theoretical[3]} & {rec.theoretical[4]} \\\\"
                )

            if not (dist_idx == len(dist_order) - 1 and c_idx == len(c_values) - 1):
                lines.append("            \\midrule")

    return lines


def render_table(body_lines: List[str], cfg: TableConfig) -> str:
    lines = [
        r"\begin{table}[htbp]",
        "",
        f"    {cfg.caption}",
        "",
        f"    {cfg.label}",
        "",
        r"    \scriptsize",
        "",
        r"    \begin{tabular}{ccccrrrrrr}",
        "",
        r"        \toprule",
        "",
        f"        {cfg.header_line1}",
        "",
        f"        {cfg.header_line2}",
        "",
        f"        {cfg.header_line3}",
        "",
        r"        \midrule",
        "",
    ]
    lines.extend(body_lines)
    lines.extend(
        [
            "",
            r"        \bottomrule",
            "",
            r"    \end{tabular}",
            "",
            r"\end{table}",
            "",
        ]
    )
    return "\n".join(lines)


def convert_one(input_xlsx: Path, output_tex: Path, cfg: TableConfig) -> None:
    rows = _read_sheet_rows(input_xlsx)
    _, records = parse_rows_to_records(rows)
    if not records:
        raise ValueError(f"No records parsed from Excel: {input_xlsx}")

    body = build_table_body(records)
    latex = render_table(body, cfg)
    output_tex.write_text(latex, encoding="utf-8")
    print(f"Wrote LaTeX table to: {output_tex}")


def main() -> None:
    default_dir = Path(__file__).resolve().parent
    parser = argparse.ArgumentParser(description="Convert Excel results to LaTeX table.")
    parser.add_argument(
        "--dataset",
        choices=sorted(TABLE_CONFIGS.keys()) + ["all"],
        default="all",
        help="Built-in dataset preset; use 'all' to generate both tables.",
    )
    parser.add_argument(
        "--xlsx",
        type=Path,
        default=None,
        help="Input .xlsx file.",
    )
    parser.add_argument(
        "--out",
        type=Path,
        default=None,
        help="Output .tex file.",
    )
    args = parser.parse_args()

    if args.dataset == "all":
        if args.xlsx is not None or args.out is not None:
            raise ValueError("--xlsx/--out cannot be used when --dataset all.")
        for key in ("bidiagonal", "spike"):
            cfg = TABLE_CONFIGS[key]
            convert_one(default_dir / cfg.xlsx, default_dir / cfg.out, cfg)
        return

    cfg = TABLE_CONFIGS[args.dataset]
    input_xlsx = args.xlsx if args.xlsx is not None else (default_dir / cfg.xlsx)
    output_tex = args.out if args.out is not None else (default_dir / cfg.out)
    convert_one(input_xlsx, output_tex, cfg)


if __name__ == "__main__":
    main()
