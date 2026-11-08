#!/usr/bin/env python3
"""qqq_run_report.py -- the run-level record of what qqq did to every sample.

Writes two files from one model:

  <prefix>.json   for machines. A stable, versioned schema (SCHEMA_VERSION, documented in
                  docs/output.md): run-level `tools` and `databases`, a `summary`, and a
                  `samples` array with per-step input/removed/output read counts and the
                  published location of the final reads and of every removal step's
                  removed reads -- or null plus a `reason` when they were not kept.
  <prefix>.html   for people. Self-contained (no network), light/dark, with stat tiles, a
                  sortable and filterable per-sample table and small inline bars.

Inputs:
  --manifest      JSON written by the workflow (workflows/qqq.nf, qqqReportManifest). It
                  carries everything only Nextflow knows: run metadata, params, database
                  paths and any checksums, the versions topic, the arguments the pipeline
                  declares, and one record per published file or tool report, already
                  resolved to its path under params.outdir.
  --reports-dir   the tool reports the manifest names, staged flat: fastp/fastplong JSON,
                  deacon/hostile/hostess JSON and bbduk stats files.

Everything numeric is parsed here from the tools' own reports rather than in Groovy, so the
parsing is covered by this module's nf-test fixtures.

Units: "reads" are sequence records. Paired-end mates count separately, which is how every
tool in the chain counts (fastp total_reads, deacon seqs_in, hostile reads_in, bbduk #Total),
so the numbers chain from step to step without halving anything.

Standard library only.
"""
import argparse
import datetime
import html
import json
import os
import platform
import re
import shlex
import sys

SCHEMA_VERSION = "1.0.0"

# Order matters: it is the order the steps run in, and the order a sample's reads flow
# through them. Each step's input should be the previous step's output.
STEP_ORDER = ["read_filter", "host_removal", "rrna_depletion", "downsampling"]
STEP_LABEL = {
    "read_filter": "Read filter",
    "host_removal": "Host removal",
    "rrna_depletion": "rRNA depletion",
    "downsampling": "Downsampling",
}

HOST_NOT_RETAINED = "host reads not retained (deacon/hostile keep none)"

# fastp/fastplong options naming a per-sample file. They are dropped from the tool-level
# argument summary so samples that ran with the same settings collapse into one line; the
# full per-run command is kept on the sample.
FASTP_IO_FLAGS = {
    "--in1", "--in2", "--out1", "--out2", "-i", "-I", "-o", "-O", "--json", "-j", "--html",
    "-h", "--failed_out", "--unpaired1", "--unpaired2", "--merged_out", "--in", "--out",
    "--report_title", "-R",
}
FASTP_IO_SWITCHES = {"--stdout"}

# deacon summary keys that are settings rather than counts, timings or file names.
DEACON_SETTINGS = ["k", "w", "abs_threshold", "rel_threshold", "prefix_length", "deplete",
                   "rename", "ordered", "check_pairs"]


# ------------------------------------------------------------------------------------------ parsing


def load_json(path):
    """Parsed JSON, or None for a missing or empty file (the module stubs write empty ones)."""
    try:
        with open(path) as fh:
            text = fh.read().strip()
    except OSError:
        return None
    if not text:
        return None
    return json.loads(text)


def parse_fastp(path):
    """fastp or fastplong JSON -> counts, filter reasons, adapter effect and the command."""
    d = load_json(path)
    if not d or "summary" not in d:
        return None
    s = d["summary"]
    before, after = s.get("before_filtering", {}), s.get("after_filtering", {})
    fr = d.get("filtering_result", {}) or {}
    ac = d.get("adapter_cutting", {}) or {}
    version = s.get("fastp_version") or s.get("fastplong_version")
    # Every *_reads key except the pass count is a reason a read was dropped. Taken
    # generically so dedup (duplicated_reads) or low-complexity filtering, if switched on in
    # config, appears without a change here.
    reasons = {k[:-len("_reads")]: v for k, v in fr.items()
               if k.endswith("_reads") and k != "passed_filter_reads"}
    return {
        "input_reads": before.get("total_reads"),
        "input_bases": before.get("total_bases"),
        "output_reads": after.get("total_reads"),
        "output_bases": after.get("total_bases"),
        "removed_by_reason": reasons,
        "adapter_trimmed_reads": ac.get("adapter_trimmed_reads"),
        "adapter_trimmed_bases": ac.get("adapter_trimmed_bases"),
        "command": d.get("command", "").strip() or None,
        "version": version,
    }


def parse_host(path, tool):
    """deacon (object, seqs_*) or hostile/hostess (one-element array, reads_*) -> counts."""
    d = load_json(path)
    if d is None:
        return None
    if isinstance(d, list):
        d = d[0] if d else {}
    if "seqs_in" in d:
        out = {
            "input_reads": d.get("seqs_in"), "output_reads": d.get("seqs_out"),
            "removed_reads": d.get("seqs_removed"), "input_bases": d.get("bp_in"),
            "output_bases": d.get("bp_out"), "removed_bases": d.get("bp_removed"),
            "settings": {k: d[k] for k in DEACON_SETTINGS if k in d},
            "index": d.get("index"), "version": d.get("version"),
        }
    else:
        out = {
            "input_reads": d.get("reads_in"), "output_reads": d.get("reads_out"),
            "removed_reads": d.get("reads_removed"), "input_bases": None,
            "output_bases": None, "removed_bases": None,
            "settings": {"aligner": d.get("aligner"), "options": d.get("options")},
            "index": d.get("index"), "version": d.get("version"),
        }
    out["tool"] = tool
    return out


def parse_bbduk(path):
    """bbduk stats= file header: #Total is records in, #Matched records split off."""
    total = matched = None
    try:
        with open(path) as fh:
            for line in fh:
                if not line.startswith("#"):
                    break
                f = line.rstrip("\n").split("\t")
                if f[0] == "#Total" and len(f) > 1:
                    total = int(f[1])
                elif f[0] == "#Matched" and len(f) > 1:
                    matched = int(f[1])
    except (OSError, ValueError):
        return None
    if total is None:
        return None
    return {"input_reads": total, "removed_reads": matched,
            "output_reads": None if matched is None else total - matched}


def fastp_arguments(command):
    """`fastp --in1 a --out1 b --thread 2 --detect_adapter_for_pe` -> `--thread 2 --detect_adapter_for_pe`."""
    if not command:
        return None
    try:
        tok = shlex.split(command)[1:]
    except ValueError:
        return command
    keep, skip = [], False
    for t in tok:
        if skip:
            skip = False
            continue
        if t in FASTP_IO_FLAGS:
            skip = True
            continue
        if t in FASTP_IO_SWITCHES:
            continue
        keep.append(t)
    return " ".join(keep)


# ------------------------------------------------------------------------------------------ model


def short_process(name):
    """QQQ:DEACON_FILTER_SMALL -> DEACON_FILTER; size tiers are one tool, not two."""
    name = name.split(":")[-1]
    return re.sub(r"_(SMALL|LARGE)(_(SHORT|LONG))?$", "", name)


def frac(num, den):
    if num is None or not den:
        return None
    return num / den


def sum_known(values):
    vals = [v for v in values if v is not None]
    return sum(vals) if vals else None


def step_skeleton(status, tool=None, reason=None):
    return {"status": status, "tool": tool, "reason": reason, "input_reads": None,
            "removed_reads": None, "output_reads": None}


def build_samples(manifest, reports_dir):
    settings = manifest.get("settings", {})
    recs_by_sample = {}
    for r in manifest.get("records", []):
        recs_by_sample.setdefault(r["sample"], []).append(r)

    def report_path(name):
        return os.path.join(reports_dir, name)

    samples = []
    tool_evidence = {"fastp": {}, "fastplong": {}, "deacon": {}, "hostile": {}, "hostess": {}}
    for sid in sorted(recs_by_sample):
        recs = recs_by_sample[sid]
        run_recs = sorted((r for r in recs if r["kind"] == "input_run"), key=lambda r: r.get("run_number") or 0)
        info = run_recs[0] if run_recs else {}
        single_end = info.get("single_end")
        warn = []
        steps = {}

        # ---- read filter (fastp / fastplong), one report per run, summed over runs
        fp_recs = sorted((r for r in recs if r["kind"] == "fastp_report"), key=lambda r: r["run"])
        if settings.get("skip_fastp"):
            steps["read_filter"] = step_skeleton("skipped", reason="skip_fastp=true")
        else:
            tool = fp_recs[0]["tool"] if fp_recs else None
            st = step_skeleton("ran", tool=tool)
            runs, parsed_all = [], True
            for r in fp_recs:
                p = parse_fastp(report_path(r["file"]))
                if p is None:
                    parsed_all = False
                    runs.append({"run": r["run"], "report": r["file"]})
                    continue
                args = fastp_arguments(p["command"])
                tool_evidence[r["tool"]].setdefault(args, set()).add(sid)
                runs.append({"run": r["run"], "report": r["file"], "command": p["command"], **{
                    k: p[k] for k in ("input_reads", "input_bases", "output_reads", "output_bases",
                                      "adapter_trimmed_reads", "adapter_trimmed_bases",
                                      "removed_by_reason")}})
            if fp_recs and parsed_all:
                for k in ("input_reads", "input_bases", "output_reads", "output_bases",
                          "adapter_trimmed_reads", "adapter_trimmed_bases"):
                    st[k] = sum_known(x.get(k) for x in runs)
                if st["input_reads"] is not None and st["output_reads"] is not None:
                    st["removed_reads"] = st["input_reads"] - st["output_reads"]
                reasons = {}
                for x in runs:
                    for k, v in (x.get("removed_by_reason") or {}).items():
                        reasons[k] = reasons.get(k, 0) + (v or 0)
                st["removed_by_reason"] = reasons
                # "Did adapter removal do anything": any read had adapter sequence cut.
                st["adapter_trimming_effective"] = (st["adapter_trimmed_reads"] or 0) > 0
            else:
                st.update({"input_bases": None, "output_bases": None, "removed_by_reason": None,
                           "adapter_trimmed_reads": None, "adapter_trimmed_bases": None,
                           "adapter_trimming_effective": None})
                if not fp_recs:
                    warn.append("read filter ran but no fastp report reached the report")
                else:
                    warn.append("fastp report empty or unparseable (stub run?); counts unknown")
            st["runs"] = runs
            steps["read_filter"] = st

        # ---- host removal
        dehoster = settings.get("dehoster")
        if dehoster == "none":
            steps["host_removal"] = step_skeleton("skipped", reason="dehoster=none")
        else:
            st = step_skeleton("ran", tool=dehoster)
            st.update({"input_bases": None, "output_bases": None, "removed_bases": None})
            hr = next((r for r in recs if r["kind"] == "host_report"), None)
            p = parse_host(report_path(hr["file"]), hr["tool"]) if hr else None
            if p:
                for k in ("input_reads", "output_reads", "removed_reads", "input_bases",
                          "output_bases", "removed_bases"):
                    st[k] = p[k]
                st["report"] = hr["file"]
                key = json.dumps({"settings": p["settings"], "index": p["index"]}, sort_keys=True)
                tool_evidence[hr["tool"]].setdefault(key, set()).add(sid)
            else:
                warn.append("host-removal report empty or missing (stub run?); counts unknown")
            steps["host_removal"] = st

        # ---- rRNA depletion
        if not settings.get("rrna_depletion"):
            steps["rrna_depletion"] = step_skeleton("skipped", reason="rrna_ref not set")
        else:
            st = step_skeleton("ran", tool="bbduk")
            rr = next((r for r in recs if r["kind"] == "rrna_report"), None)
            p = parse_bbduk(report_path(rr["file"])) if rr else None
            if p:
                st.update(p)
                st["report"] = rr["file"]
                if single_end is False:
                    # #Matched counts reads whose own k-mers hit the reference, but bbduk
                    # splits pairs together: one matching mate removes both. So for paired-end
                    # data the records actually removed run slightly higher (216 against 210
                    # on the test sample; 0.27-0.28% on the validation library). The exact
                    # figure is only in bbduk's stderr, which the module does not keep.
                    st["exact"] = False
                    st["note"] = ("paired-end: from bbduk #Matched, which counts per-read matches; bbduk "
                                  "removes the whole pair, so removed_reads is a slight undercount")
            else:
                warn.append("bbduk stats file missing or unparseable; counts unknown")
            steps["rrna_depletion"] = st

        # ---- downsampling (rasusa). It writes no report, so the count is derived: `rasusa
        # reads --num` measures its input first and keeps exactly the cap, or everything
        # when the input is under it. The module's nf-test pins that contract.
        max_reads = settings.get("max_reads")
        prev_out = None
        for s in ("rrna_depletion", "host_removal", "read_filter"):
            if steps[s]["status"] == "ran":
                prev_out = steps[s]["output_reads"]
                break
        if not max_reads:
            steps["downsampling"] = step_skeleton("skipped", reason="max_reads not set")
        else:
            went = any(r["kind"] == "reads" and r["step"] == "downsampling" for r in recs)
            cap = max_reads * (1 if single_end else 2)
            if went:
                st = step_skeleton("ran", tool="rasusa")
                st["input_reads"] = prev_out
                st["output_reads"] = None if prev_out is None else min(prev_out, cap)
                st["removed_reads"] = None if prev_out is None else prev_out - st["output_reads"]
                st["derived"] = True
                st["cap_reads"] = cap
            else:
                st = step_skeleton("bypassed", tool="rasusa",
                                   reason="at or under max_reads; not routed to rasusa")
                st["input_reads"] = prev_out
                st["output_reads"] = prev_out
                st["removed_reads"] = 0 if prev_out is not None else None
                st["cap_reads"] = cap
            steps["downsampling"] = st

        # ---- chain checks: each step should start where the previous one ended
        last_out, last_name = None, None
        for s in STEP_ORDER:
            st = steps[s]
            if st["status"] not in ("ran", "bypassed"):
                continue
            if last_out is not None and st["input_reads"] is not None and st["input_reads"] != last_out:
                warn.append(f"{s} input ({st['input_reads']:,}) differs from {last_name} output ({last_out:,})")
            if st["output_reads"] is not None:
                last_out, last_name = st["output_reads"], s

        input_reads = input_bases = None
        for s in STEP_ORDER:
            if steps[s]["status"] == "ran":
                input_reads = steps[s]["input_reads"]
                input_bases = steps[s].get("input_bases")
                break
        output_reads = None
        for s in reversed(STEP_ORDER):
            if steps[s]["status"] in ("ran", "bypassed"):
                output_reads = steps[s]["output_reads"]
                break
        zero = None if output_reads is None else output_reads == 0
        if zero:
            warn.append("no reads left after processing")

        samples.append({
            "id": sid,
            "single_end": single_end,
            "platform": info.get("platform"),
            "read_type": info.get("read_type"),
            "runs": [{"id": r["run"], "input_files": r.get("files", [])} for r in run_recs],
            "input_reads": input_reads,
            "input_bases": input_bases,
            "output_reads": output_reads,
            "retained_fraction": frac(output_reads, input_reads),
            "zero_reads": zero,
            "steps": steps,
            "outputs": build_outputs(recs, steps, settings),
            "warnings": warn,
        })
    return samples, tool_evidence


def files_of(recs, kind, **match):
    out = []
    for r in recs:
        if r["kind"] == kind and all(r.get(k) == v for k, v in match.items()):
            out.extend(r.get("files", []))
    return out


def single_end_of(recs):
    return next((r.get("single_end") for r in recs if r["kind"] == "input_run"), None)


def located(files, reason=None, **extra):
    if files:
        return {"files": files, "reason": None, **extra}
    return {"files": None, "reason": reason, **extra}


def build_outputs(recs, steps, settings):
    # Final reads: the highest-priority producer the workflow saw for this sample. The
    # priorities follow the order the steps run in, so a sample under the max_reads cap
    # (absent from rasusa's output) falls back to whatever it last passed through.
    producers = [r for r in recs if r["kind"] == "reads"]
    clean = None
    if producers:
        best = max(producers, key=lambda r: r["priority"])
        if best.get("published"):
            clean = {"step": best["step"], "files": best["files"], "reason": None}
        else:
            clean = {"step": best["step"], "files": None,
                     "reason": "no step rewrote these reads, so nothing was published; they are the samplesheet inputs",
                     "source_files": best.get("source_files")}
    removed = {}

    st = steps["read_filter"]
    if st["status"] != "ran":
        removed["read_filter"] = located(None, f"step skipped ({st['reason']})")
    elif not settings.get("save_failed_reads"):
        removed["read_filter"] = located(None, "save_failed_reads=false; fastp wrote no removed reads")
    else:
        removed["read_filter"] = located(files_of(recs, "fastp_failed"),
                                         "save_failed_reads=true but no failed-read files were published")
        if st["tool"] == "fastp" and removed["read_filter"]["files"] and not single_end_of(recs):
            # Measured on fastp 1.3.6: a pair in which BOTH mates fail is counted as removed
            # but written to neither --failed_out nor --unpaired1/2. On the test sample 32 of
            # 373 removed pairs were missing. Single-end fastp and fastplong write every read.
            removed["read_filter"]["complete"] = False
            removed["read_filter"]["note"] = (
                "paired-end fastp does not write pairs in which both mates fail, so these files can "
                "hold fewer reads than removed_reads")

    st = steps["host_removal"]
    if st["status"] != "ran":
        removed["host_removal"] = located(None, f"step skipped ({st['reason']})")
    elif st["tool"] == "hostess":
        removed["host_removal"] = located(files_of(recs, "host_bam"),
                                          "host_bam=false; hostess was not asked to keep removed reads"
                                          if not settings.get("host_bam") else
                                          "host_bam=true but no host BAM was published for this sample")
    else:
        removed["host_removal"] = located(None, HOST_NOT_RETAINED)

    st = steps["rrna_depletion"]
    if st["status"] != "ran":
        removed["rrna_depletion"] = located(None, f"step skipped ({st['reason']})")
    else:
        removed["rrna_depletion"] = located(files_of(recs, "rrna_matched"),
                                            "bbduk wrote no matched reads for this sample")

    st = steps["downsampling"]
    if st["status"] == "skipped":
        removed["downsampling"] = located(None, f"step skipped ({st['reason']})")
    elif st["status"] == "bypassed":
        removed["downsampling"] = located(None, "nothing removed; sample was under max_reads")
    else:
        removed["downsampling"] = located(
            None, "rasusa writes only the reads it keeps; the complete pre-downsampling reads "
                  "are the previous step's published output")
    return {"clean_reads": clean, "removed_reads": removed}


def build_tools(manifest, evidence):
    by_tool = {}
    for v in manifest.get("versions", []):
        t = by_tool.setdefault(v["tool"], {"name": v["tool"], "versions": set(), "processes": set()})
        t["versions"].add(str(v["version"]) if v["version"] not in (None, "") else "")
        t["processes"].add(short_process(v["process"]))
    declared = {}
    for d in manifest.get("declared_arguments", []):
        declared.setdefault(d["tool"], []).append(d)
    tools = []
    for name in sorted(by_tool, key=str.lower):
        t = by_tool[name]
        versions = sorted(x for x in t["versions"] if x)
        entry = {"name": name, "version": versions[0] if len(versions) == 1 else (versions or [None])[0],
                 "versions": versions, "processes": sorted(t["processes"]),
                 "arguments": [], "arguments_source": None}
        if not versions:
            entry["version"] = None
            entry["version_note"] = "the module's version command returned nothing"
        ev = evidence.get(name)
        if name in ("fastp", "fastplong") and ev:
            entry["arguments"] = [{"arguments": a, "samples": len(s)} for a, s in sorted(ev.items(), key=lambda kv: str(kv[0]))]
            entry["arguments_source"] = "tool_report"
        elif name in ("deacon", "hostile", "hostess") and ev:
            entry["arguments"] = [{**json.loads(k), "samples": len(s)} for k, s in sorted(ev.items())]
            entry["arguments_source"] = "tool_report"
        if name in declared:
            for d in declared[name]:
                entry["arguments"].append({"process": d["process"], "arguments": d["arguments"]})
            entry["arguments_source"] = ("tool_report+pipeline_config" if entry["arguments_source"]
                                         else "pipeline_config")
        tools.append(entry)
    tools.append({"name": "python", "version": platform.python_version(),
                  "versions": [platform.python_version()], "processes": ["QQQ_RUN_REPORT"],
                  "arguments": [], "arguments_source": None,
                  "version_note": "recorded by the report itself; QQQ_RUN_REPORT reads the versions "
                                  "topic, so it cannot also emit to it"})
    return tools


def build_databases(manifest, evidence):
    dbs = []
    for d in manifest.get("databases", []):
        d = dict(d)
        tool = d.get("tool")
        # The per-sample reports record what the tool saw at runtime: deacon its k, w and
        # thresholds, hostile/hostess the aligner and index name. Attach the distinct values.
        if tool in evidence and evidence[tool] and d.get("role", "").startswith("host"):
            seen = [json.loads(k) for k in sorted(evidence[tool])]
            d.setdefault("parameters", {})
            d["parameters"]["reported_by_tool"] = seen if len(seen) > 1 else seen[0]
        dbs.append(d)
    return dbs


def build_summary(samples):
    def tot(key, step=None):
        vals = [(s["steps"][step][key] if step else s[key]) for s in samples
                if (s["steps"][step]["status"] in ("ran", "bypassed") if step else True)]
        return sum_known(vals)
    inp, out = tot("input_reads"), tot("output_reads")
    removed = {s: tot("removed_reads", s) for s in STEP_ORDER}
    trimmed = sum_known(s["steps"]["read_filter"].get("adapter_trimmed_reads")
                        for s in samples if s["steps"]["read_filter"]["status"] == "ran")
    fp_in = tot("input_reads", "read_filter")
    reasons = {}
    for s in samples:
        for k, v in (s["steps"]["read_filter"].get("removed_by_reason") or {}).items():
            reasons[k] = reasons.get(k, 0) + (v or 0)
    return {
        "samples": len(samples),
        "samples_zero_reads": sum(1 for s in samples if s["zero_reads"]),
        "samples_zero_reads_ids": [s["id"] for s in samples if s["zero_reads"]],
        "samples_with_warnings": sum(1 for s in samples if s["warnings"]),
        "input_reads": inp,
        "input_bases": tot("input_bases"),
        "output_reads": out,
        "retained_fraction": frac(out, inp),
        "removed_reads": removed,
        "removed_fraction": {k: frac(v, inp) for k, v in removed.items()},
        "read_filter_removed_by_reason": reasons,
        "adapter_trimmed_reads": trimmed,
        "adapter_trimmed_fraction": frac(trimmed, fp_in),
        "adapter_trimmed_bases": sum_known(s["steps"]["read_filter"].get("adapter_trimmed_bases")
                                           for s in samples if s["steps"]["read_filter"]["status"] == "ran"),
    }


def build(manifest, reports_dir):
    samples, evidence = build_samples(manifest, reports_dir)
    return {
        "schema_version": SCHEMA_VERSION,
        "generated_at": datetime.datetime.now(datetime.timezone.utc).isoformat(timespec="seconds"),
        "generator": {"script": "qqq_run_report.py", "python": platform.python_version()},
        "units": "reads are sequence records; paired-end mates are counted separately",
        "step_order": STEP_ORDER,
        "run": manifest.get("run", {}),
        "tools": build_tools(manifest, evidence),
        "databases": build_databases(manifest, evidence),
        "summary": build_summary(samples),
        "outputs": manifest.get("run_outputs", {}),
        "samples": samples,
    }


# ------------------------------------------------------------------------------------------ html

def esc(s):
    return html.escape("" if s is None else str(s), quote=True)


def human(n):
    if n is None:
        return "–"
    n = float(n)
    for div, suf in ((1e9, "G"), (1e6, "M"), (1e3, "k")):
        if abs(n) >= div:
            return f"{n / div:.1f}{suf}"
    return f"{int(n)}"


def pct(x, digits=None):
    if x is None:
        return "–"
    if digits is None:
        digits = 2 if 0 < abs(x) < 0.01 else 1
    return f"{x * 100:.{digits}f}%"


def kv(d):
    """{'k': 31, 'settings': {'w': 15}} -> 'k=31 w=15'; nested settings flattened, lists comma-joined."""
    out = []
    for k, v in sorted(d.items()):
        if isinstance(v, dict):
            out.append(kv(v))
        elif isinstance(v, list) and v and all(isinstance(x, dict) for x in v):
            out.append(" | ".join(kv(x) for x in v))
        elif isinstance(v, list):
            out.append(f"{k}={','.join(map(str, v)) or '(none)'}")
        else:
            out.append(f"{k}={json.dumps(v) if isinstance(v, bool) or v is None else v}")
    return " ".join(x for x in out if x)


SEGMENTS = [("read_filter", "seg-filter"), ("host_removal", "seg-host"),
            ("rrna_depletion", "seg-rrna"), ("downsampling", "seg-down")]


def comp_bar(s, width=150):
    """Stacked bar of where a sample's input reads went: removed per step, then retained."""
    inp = s["input_reads"]
    if not inp:
        return '<span class="sub">–</span>'
    x, parts = 0.0, []
    for step, cls in SEGMENTS:
        v = s["steps"][step]["removed_reads"] or 0
        w = width * v / inp
        if w > 0:
            parts.append(f'<rect x="{x:.1f}" width="{w:.1f}" height="10" class="{cls}"/>')
        x += w
    keep = width * (s["output_reads"] or 0) / inp
    parts.append(f'<rect x="{x:.1f}" width="{keep:.1f}" height="10" class="seg-keep"/>')
    title = "; ".join(f"{STEP_LABEL[st]} {pct(frac(s['steps'][st]['removed_reads'], inp))}"
                      for st, _ in SEGMENTS if s["steps"][st]["removed_reads"])
    return (f'<svg width="{width}" height="10" class="bar" role="img" aria-label="{esc(title)}">'
            f'<title>{esc(title or "nothing removed")}</title>'
            f'<rect width="{width}" height="10" rx="2" class="track"/>{"".join(parts)}</svg>')


def count_cell(n, of, title=""):
    """A count with its share of the input underneath; sorts on the share when there is one."""
    share = frac(n, of)
    key = share if share is not None else (n if n is not None else -1)
    sub = f'<div class="sub">{pct(share)}</div>' if share is not None else ""
    return (f'<td class="num" data-v="{key}" title="{esc(title or (f"{n:,}" if n is not None else ""))}">'
            f"{human(n)}{sub}</td>")


def step_cell(s, step):
    st = s["steps"][step]
    if st["status"] == "skipped":
        return '<td class="num na" data-v="-1">off</td>'
    if st["status"] == "bypassed":
        return '<td class="num na" data-v="0" title="under max_reads">–</td>'
    title = f"{st['removed_reads']:,} of {st['input_reads']:,} reads" if st["removed_reads"] is not None and st["input_reads"] else ""
    if step == "read_filter" and st.get("removed_by_reason"):
        title += " · " + ", ".join(f"{k}: {v:,}" for k, v in st["removed_by_reason"].items() if v)
    if st.get("note"):
        title += " · " + st["note"]
    return count_cell(st["removed_reads"], s["input_reads"], title)


def path_list(loc):
    if loc is None:
        return '<span class="sub">–</span>'
    if loc.get("files"):
        note = f'<div class="reason">{esc(loc["note"])}</div>' if loc.get("note") else ""
        return "".join(f'<div class="path" title="{esc(f["path"])}">{esc(f["relative_path"])}</div>'
                       for f in loc["files"]) + note
    return f'<span class="reason">{esc(loc.get("reason"))}</span>'


CSS = """
:root{--bg:#fbfbfa;--fg:#1d1f21;--mut:#6b7075;--line:#e3e4e2;--track:#ecedeb;--card:#fff;--code:#f2f2f0;
--filter:#8a8f94;--host:#c47a2c;--rrna:#7b68b5;--down:#b8c4cf;--keep:#3d8b5a;--warn:#b4413d;--accent:#3f6e9e}
@media (prefers-color-scheme:dark){:root:not([data-theme="light"]){--bg:#16181a;--fg:#e6e6e3;--mut:#9aa0a6;--line:#2c3033;
--track:#2a2e31;--card:#1d2023;--code:#22262a;--filter:#80868c;--host:#d9934a;--rrna:#9a8ad0;--down:#5d6b78;--keep:#5aae78;
--warn:#e0706b;--accent:#7fa3cc;color-scheme:dark}}
:root[data-theme="dark"]{--bg:#16181a;--fg:#e6e6e3;--mut:#9aa0a6;--line:#2c3033;--track:#2a2e31;--card:#1d2023;
--code:#22262a;--filter:#80868c;--host:#d9934a;--rrna:#9a8ad0;--down:#5d6b78;--keep:#5aae78;--warn:#e0706b;
--accent:#7fa3cc;color-scheme:dark}
*{box-sizing:border-box}
body{margin:0;background:var(--bg);color:var(--fg);font:13px/1.45 system-ui,-apple-system,Segoe UI,sans-serif;
padding:24px clamp(16px,3vw,40px)}
h1{font-size:20px;margin:0 0 4px}
h2{font-size:15px;margin:30px 0 8px;border-bottom:1px solid var(--line);padding-bottom:4px}
.sub{color:var(--mut);font-size:11px}
.head{display:flex;justify-content:space-between;gap:16px;align-items:flex-start;flex-wrap:wrap}
button.theme{background:var(--card);color:var(--fg);border:1px solid var(--line);border-radius:4px;padding:3px 10px;cursor:pointer;font:inherit}
.tiles{display:grid;grid-template-columns:repeat(auto-fill,minmax(150px,1fr));gap:10px;margin:14px 0}
.tile{background:var(--card);border:1px solid var(--line);border-radius:6px;padding:8px 12px;border-top:3px solid var(--c,var(--line))}
.tile b{display:block;font-size:20px;font-variant-numeric:tabular-nums}
.tile .sub{display:block}
.tw{overflow-x:auto;max-width:100%}
table{border-collapse:collapse;font-variant-numeric:tabular-nums}
th{text-align:left;font-weight:600;color:var(--mut);font-size:11px;padding:4px 8px;border-bottom:1px solid var(--line);
position:sticky;top:0;background:var(--bg);white-space:nowrap}
#t th{cursor:pointer;user-select:none} #t th[aria-sort=ascending]::after{content:" ▲"} #t th[aria-sort=descending]::after{content:" ▼"}
td{padding:4px 8px;vertical-align:top;border-bottom:1px solid var(--line)}
td.num{text-align:right;white-space:nowrap} td.na{color:var(--mut)}
td.id{font-weight:600;white-space:nowrap}
.bar .track{fill:var(--track)} .seg-filter{fill:var(--filter)} .seg-host{fill:var(--host)} .seg-rrna{fill:var(--rrna)}
.seg-down{fill:var(--down)} .seg-keep{fill:var(--keep)}
.legend{display:flex;flex-wrap:wrap;gap:12px;align-items:center;font-size:11px;color:var(--mut)}
.sw{display:inline-block;width:10px;height:10px;border-radius:2px;margin-right:4px;vertical-align:-1px}
.flag{font-size:10px;font-weight:600;padding:1px 6px;border-radius:9px;color:#fff;background:var(--warn);white-space:nowrap}
.controls{display:flex;flex-wrap:wrap;gap:14px;align-items:center;margin:8px 0}
.controls input[type=search]{padding:4px 8px;border:1px solid var(--line);border-radius:4px;background:var(--card);color:var(--fg);min-width:220px}
code,.path{font:11px/1.4 ui-monospace,SFMono-Regular,Menlo,Consolas,monospace;word-break:break-all}
code{background:var(--code);padding:1px 4px;border-radius:3px}
.reason{color:var(--mut);font-size:11px;font-style:italic}
.kv td:first-child{color:var(--mut);white-space:nowrap}
.note{color:var(--mut);font-size:12px;max-width:900px}
ul.w{margin:0;padding-left:16px;color:var(--warn);font-size:11px}
"""

JS = """
(function(){
const root=document.documentElement,btn=document.getElementById('theme');
function cur(){return root.dataset.theme||(matchMedia('(prefers-color-scheme: dark)').matches?'dark':'light');}
btn.addEventListener('click',()=>{root.dataset.theme=cur()==='dark'?'light':'dark';});
const t=document.getElementById('t');if(!t)return;
const body=t.tBodies[0],rows=[...body.rows],q=document.getElementById('q'),z=document.getElementById('zero'),
w=document.getElementById('warn'),n=document.getElementById('nshown');
function apply(){const s=q.value.toLowerCase();let k=0;rows.forEach(r=>{let show=true;
if(s&&!r.dataset.id.includes(s))show=false;if(z.checked&&r.dataset.zero!=='1')show=false;
if(w.checked&&r.dataset.warn!=='1')show=false;r.hidden=!show;if(show)k++;});n.textContent=k+' of '+rows.length+' samples';}
[q,z,w].forEach(e=>e.addEventListener('input',apply));apply();
[...t.tHead.rows[0].cells].forEach((th,i)=>th.addEventListener('click',()=>{
const asc=th.getAttribute('aria-sort')!=='ascending';[...t.tHead.rows[0].cells].forEach(c=>c.removeAttribute('aria-sort'));
th.setAttribute('aria-sort',asc?'ascending':'descending');const num=th.dataset.type==='num';
rows.sort((a,b)=>{const x=a.cells[i],y=b.cells[i];let r;
if(num){r=parseFloat(x.dataset.v??-1)-parseFloat(y.dataset.v??-1);}else{r=x.textContent.localeCompare(y.textContent);}
return asc?r:-r;});rows.forEach(r=>body.appendChild(r));}));
})();
"""


def render(rep):
    run, summ, samples = rep["run"], rep["summary"], rep["samples"]
    title = "qqq run report"
    P = []
    commit = run.get("commit") or "unknown commit"
    P.append('<div class="head"><div>')
    P.append(f"<h1>{esc(title)} · {esc(run.get('run_name', ''))}</h1>")
    P.append(f'<div class="sub">{esc(run.get("pipeline", ""))} {esc(run.get("pipeline_version", ""))} · '
             f'{esc(commit[:12])}{" (dirty)" if run.get("commit_dirty") else ""} · outdir <code>{esc(run.get("outdir", ""))}</code> · '
             f'generated {esc(rep["generated_at"])} · schema {esc(rep["schema_version"])}</div>')
    P.append('</div><button class="theme" id="theme" type="button">Light / dark</button></div>')

    inp = summ["input_reads"]
    rm = summ["removed_reads"]
    tiles = [
        (summ["samples"], "samples", None, "var(--accent)"),
        (human(inp), "input reads", f"{human(summ['input_bases'])} bases" if summ["input_bases"] else None, "var(--accent)"),
        (pct(summ["adapter_trimmed_fraction"]), "reads adapter-trimmed",
         f"{human(summ['adapter_trimmed_reads'])} reads · {human(summ['adapter_trimmed_bases'])} bases"
         if summ["adapter_trimmed_reads"] is not None else "read filter off", "var(--line)"),
        (human(rm["read_filter"]), "removed by read filter", pct(summ["removed_fraction"]["read_filter"]), "var(--filter)"),
        (human(rm["host_removal"]), "removed as host", pct(summ["removed_fraction"]["host_removal"]), "var(--host)"),
        (human(rm["rrna_depletion"]), "removed as rRNA", pct(summ["removed_fraction"]["rrna_depletion"]), "var(--rrna)"),
        (human(rm["downsampling"]), "removed by downsampling", pct(summ["removed_fraction"]["downsampling"]), "var(--down)"),
        (human(summ["output_reads"]), "final reads", f"{pct(summ['retained_fraction'])} of input", "var(--keep)"),
        (summ["samples_zero_reads"], "samples with zero reads left",
         ", ".join(summ["samples_zero_reads_ids"][:5]) or None, "var(--warn)" if summ["samples_zero_reads"] else "var(--line)"),
    ]
    P.append('<div class="tiles">' + "".join(
        f'<div class="tile" style="--c:{c}"><b>{esc(v)}</b>{esc(lbl)}'
        + (f'<span class="sub">{esc(sub)}</span>' if sub else "") + "</div>"
        for v, lbl, sub, c in tiles) + "</div>")
    reasons = summ["read_filter_removed_by_reason"]
    if reasons:
        P.append('<p class="note">Read filter removals by reason: ' + ", ".join(
            f"{esc(k.replace('_', ' '))} {v:,}" for k, v in sorted(reasons.items(), key=lambda kv: -kv[1])) + ".</p>")

    # ---- per-sample table
    P.append("<h2>Samples</h2>")
    P.append('<div class="controls"><input type="search" id="q" placeholder="filter sample…" aria-label="Filter samples">'
             '<label><input type="checkbox" id="zero"> zero reads left</label>'
             '<label><input type="checkbox" id="warn"> with warnings</label><span class="sub" id="nshown"></span>'
             '<div class="legend">'
             + "".join(f'<span><span class="sw" style="background:var(--{c})"></span>{lbl}</span>'
                       for c, lbl in (("filter", "read filter"), ("host", "host"), ("rrna", "rRNA"),
                                      ("down", "downsampling"), ("keep", "retained")))
             + "</div></div>")
    P.append('<p class="note">Counts are reads (paired-end mates counted separately); the small figure is the '
             'share of the sample\'s input. Click a column to sort. Hover a count for the exact number.</p>')
    P.append('<div class="tw"><table id="t"><thead><tr><th>Sample</th><th data-type="num">Input reads</th>'
             '<th data-type="num">Adapter-trimmed</th><th data-type="num">Read filter</th><th data-type="num">Host</th>'
             '<th data-type="num">rRNA</th><th data-type="num">Downsampled</th><th data-type="num">Final reads</th>'
             '<th>Where the input went</th><th>Flags</th></tr></thead><tbody>')
    for s in samples:
        rf = s["steps"]["read_filter"]
        at = rf.get("adapter_trimmed_reads") if rf["status"] == "ran" else None
        at_cell = (count_cell(at, rf["input_reads"], f"{at:,} reads, {rf.get('adapter_trimmed_bases') or 0:,} bases trimmed")
                   if at is not None else '<td class="num na" data-v="-1">–</td>')
        flags = ""
        if s["zero_reads"]:
            flags += '<span class="flag">zero reads</span> '
        if s["warnings"]:
            flags += '<ul class="w">' + "".join(f"<li>{esc(w)}</li>" for w in s["warnings"]) + "</ul>"
        P.append(
            f'<tr data-id="{esc(s["id"].lower())}" data-zero="{1 if s["zero_reads"] else 0}" '
            f'data-warn="{1 if s["warnings"] else 0}">'
            f'<td class="id">{esc(s["id"])}<div class="sub">{"SE" if s["single_end"] else "PE"} · '
            f'{esc(s.get("platform") or "")} · {len(s["runs"])} run{"s" if len(s["runs"]) != 1 else ""}</div></td>'
            + count_cell(s["input_reads"], None, f"{s['input_reads']:,} reads" if s["input_reads"] is not None else "")
            + at_cell + step_cell(s, "read_filter") + step_cell(s, "host_removal")
            + step_cell(s, "rrna_depletion") + step_cell(s, "downsampling")
            + count_cell(s["output_reads"], s["input_reads"])
            + f"<td>{comp_bar(s)}</td><td>{flags}</td></tr>")
    P.append("</tbody></table></div>")

    # ---- processing
    P.append("<h2>What processed this run</h2>")
    P.append('<p class="note">Versions come from the pipeline\'s versions topic. Arguments marked '
             '<i>tool report</i> were recorded by the tool itself at runtime; <i>pipeline config</i> means the '
             "pipeline's own conf/modules.config or a module input, so a launch-time <code>-c</code> override of "
             "<code>ext.args</code> would not show there.</p>")
    P.append('<div class="tw"><table><thead><tr><th>Tool</th><th>Version</th><th>Process</th><th>Arguments</th>'
             "<th>Source</th></tr></thead><tbody>")
    for t in rep["tools"]:
        args = []
        for a in t["arguments"]:
            if "arguments" in a and len(a) <= 3:
                lab = f' <span class="sub">({a["samples"]} samples)</span>' if "samples" in a else ""
                args.append(f'<div><code>{esc(a["arguments"]) or "(none)"}</code>{lab}</div>')
            else:
                shown = {k: v for k, v in a.items() if k != "samples"}
                args.append(f'<div><code>{esc(kv(shown))}</code>'
                            f'<span class="sub"> ({a.get("samples", "?")} samples)</span></div>')
        ver = esc(", ".join(t["versions"]) or "–")
        if t.get("version_note"):
            ver += f'<div class="sub">{esc(t["version_note"])}</div>'
        src = (t["arguments_source"] or "–").replace("_", " ").replace("+", " + ")
        P.append(f"<tr><td><b>{esc(t['name'])}</b></td><td>{ver}</td><td>{esc(', '.join(t['processes']))}</td>"
                 f"<td>{''.join(args) or '<span class=sub>–</span>'}</td><td class=\"sub\">{esc(src)}</td></tr>")
    P.append("</tbody></table></div>")

    P.append("<h2>Databases and references</h2>")
    P.append('<div class="tw"><table><thead><tr><th>Role</th><th>Name</th><th>Path</th><th>Checksums / size</th>'
             "<th>Parameters</th></tr></thead><tbody>")
    if not rep["databases"]:
        P.append('<tr><td colspan="5" class="sub">No database or reference was used.</td></tr>')
    for d in rep["databases"]:
        ck = []
        for k in ("md5", "sha256"):
            if d.get(k):
                ck.append(f"<div>{k} <code>{esc(d[k])}</code></div>")
        if d.get("checksum_source"):
            ck.append(f'<div class="sub">from {esc(d["checksum_source"])}</div>')
        if d.get("size_bytes") is not None:
            ck.append(f'<div class="sub">{d["size_bytes"]:,} bytes</div>')
        if not ck:
            ck.append(f'<span class="reason">{esc(d.get("checksum_note") or "none available")}</span>')
        params = d.get("parameters") or {}
        P.append(f"<tr><td>{esc(d.get('role'))}<div class=\"sub\">{esc(d.get('tool'))}</div></td>"
                 f"<td><b>{esc(d.get('name'))}</b></td><td class=\"path\">{esc(d.get('path'))}</td>"
                 f"<td>{''.join(ck)}</td><td><code>{esc(kv(params)) if params else '–'}</code></td></tr>")
    P.append("</tbody></table></div>")

    # ---- outputs
    P.append("<h2>Where the reads are</h2>")
    P.append('<p class="note">Final cleaned reads and, for each removal step, where its removed reads were saved. '
             f'Paths are relative to the outdir, <code>{esc(run.get("outdir", ""))}</code>; hover one for the full path, '
             'and the JSON carries both. Host reads are deliberately not kept unless hostess was run with '
             '<code>--host_bam</code>.</p>')
    P.append('<div class="tw"><table><thead><tr><th>Sample</th><th>Final reads</th>'
             + "".join(f"<th>Removed: {STEP_LABEL[s]}</th>" for s in STEP_ORDER) + "</tr></thead><tbody>")
    for s in samples:
        o = s["outputs"]
        cr = o["clean_reads"]
        cr_html = path_list(cr) + (f'<div class="sub">from {esc(STEP_LABEL.get(cr["step"], cr["step"]))}</div>'
                                   if cr and cr.get("step") in STEP_LABEL else "")
        P.append(f'<tr><td class="id">{esc(s["id"])}</td><td>{cr_html}</td>'
                 + "".join(f"<td>{path_list(o['removed_reads'][st])}</td>" for st in STEP_ORDER) + "</tr>")
    P.append("</tbody></table></div>")

    # ---- run
    P.append("<h2>Run</h2>")
    keys = ["pipeline", "pipeline_version", "commit", "commit_source", "revision", "repository", "run_name",
            "session_id", "started", "nextflow_version", "profile", "container_engine", "outdir", "launch_dir",
            "work_dir", "project_dir", "command_line"]
    P.append('<div class="tw"><table class="kv"><tbody>')
    for k in keys:
        if run.get(k) not in (None, ""):
            P.append(f"<tr><td>{esc(k)}</td><td><code>{esc(run[k])}</code></td></tr>")
    for k, v in (rep.get("outputs") or {}).items():
        if v:
            P.append(f"<tr><td>{esc(k)}</td><td><code>{esc(v.get('path'))}</code></td></tr>")
    P.append("</tbody></table></div>")
    if run.get("params"):
        P.append("<h2>Parameters</h2>")
        P.append('<div class="tw"><table class="kv"><tbody>' + "".join(
            f"<tr><td>{esc(k)}</td><td><code>{esc(json.dumps(v))}</code></td></tr>"
            for k, v in sorted(run["params"].items())) + "</tbody></table></div>")

    body = "\n".join(P)
    return ('<!DOCTYPE html><html lang="en"><head><meta charset="utf-8">'
            '<meta name="viewport" content="width=device-width,initial-scale=1">'
            f"<title>qqq run report</title><style>{CSS}</style></head><body>\n{body}\n"
            f"<script>{JS}</script></body></html>\n")


# ------------------------------------------------------------------------------------------ main


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--manifest", required=True)
    ap.add_argument("--reports-dir", default="reports")
    ap.add_argument("--prefix", default="qqq_run_report")
    a = ap.parse_args()

    with open(a.manifest) as fh:
        manifest = json.load(fh)
    rep = build(manifest, a.reports_dir)
    with open(f"{a.prefix}.json", "w") as fh:
        json.dump(rep, fh, indent=2, sort_keys=False)
        fh.write("\n")
    with open(f"{a.prefix}.html", "w") as fh:
        fh.write(render(rep))
    s = rep["summary"]
    print(f"{a.prefix}: {s['samples']} samples, {s['samples_zero_reads']} with zero reads left, "
          f"{s['samples_with_warnings']} with warnings", file=sys.stderr)


if __name__ == "__main__":
    main()
