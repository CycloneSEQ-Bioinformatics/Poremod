# poremod validate

MM/ML structural validation: checks the modification-tag structure of a modbam read by read, and prints a graded report grouped into four layers, to help spot corrupted or malformed modbams before analysis.

Modification codes (including the `mm_code_unknown` and `mm_code_base_mismatch` checks) are described in the [Modification codes](../../README.md#modification-codes) section of the README.

## Usage

```bash
poremod validate --input in.bam [--json] [--all]
```

## Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `--input PATH` | yes | Input modbam path (no `.bai` index needed) |
| `--json` | no | Emit the full layered report as JSON |
| `--all` | no | Scan the entire file; default is the first 10000 records |

**BREAKING:** files with more than 10000 records are no longer fully scanned unless `--all` is set.

## Checks

Layers are keyed by the object each one validates: which modification tags exist (`tag_set`), whether the MM string parses (`mm_grammar`), whether the parse result agrees with SEQ (`mm_semantics`), and the ML probability array itself (`ml`). Every issue is reported as a `(check_id, severity, anomaly_type)` triple. A check's `severity` is a fixed property of the check, not of an individual finding.

### `tag_set` — which modification tags exist

| check_id | severity | anomaly_type | Fires when |
|----------|----------|--------------|------------|
| `mm_missing` | error | `wrong_type` | MM is present but is not a `Z` (string) aux field |
| `mm_missing` | error | `ml_without_mm` | ML is present and neither `MM` nor `Mm` is |
| `ml_missing` | error | `mm_without_ml` | A usable MM string is present and neither `ML` nor `Ml` is |
| `tag_duplicate` | error | `duplicate_mm` | The read carries more than one MM tag |
| `tag_duplicate` | error | `duplicate_ml` | The read carries more than one ML tag |
| `tag_alias_conflict` | error | `alias_mm` | The read carries both `MM` and `Mm` |
| `tag_alias_conflict` | error | `alias_ml` | The read carries both `ML` and `Ml` |

`MM`/`Mm` and `ML`/`Ml` are the two spellings SAMtags allows. Reading is first-match-wins, in `MM` before `Mm` and `ML` before `Ml` order, so one lowercase alias on its own is legal and reports nothing.

### `mm_grammar` — whether the MM string parses

The raw MM string is scanned against the SAMtags §1.7 grammar `[ACGTUN][-+]([a-z]+|[0-9]+)[.?]?(,[0-9]+)*;`. Every violation in the string is collected, not just the first. A reason that quotes the offending text abbreviates it the same way the `ml` layer does: text longer than 40 characters is cut at a comma and closed with `...`.

| check_id | severity | anomaly_type | Fires when |
|----------|----------|--------------|------------|
| `mm_base_invalid` | error | `invalid_base` | A segment is shorter than two characters, or its first character is not one of `ACGTUN` |
| `mm_strand_invalid` | error | `invalid_strand` | A segment's second character is not `+` or `-` |
| `mm_code_invalid` | error | `invalid_code` | A segment's modification code is empty, mixes letters with digits, or is a ChEBI code too large for `u32` |
| `mm_marker_invalid` | error | `invalid_marker` | `?` or `.` appears anywhere but at the end of a segment's code part |
| `mm_no_terminator` | error | `missing_terminator` | The MM string does not end with `;` |
| `mm_empty_segment` | error | `empty_segment` | The MM string contains an empty segment (`;;`) |
| `mm_delta_invalid` | error | `invalid_delta` | A comma-separated delta is not a non-negative integer |

### `mm_semantics` — whether the parse result agrees with the record

| check_id | severity | anomaly_type | Fires when |
|----------|----------|--------------|------------|
| `mm_position_unresolved` | error | `position_unresolved` | A segment's delta list resolves past the end of SEQ |
| `mm_code_unknown` | warning | `unsupported` | The mod code is not in the known-code table |
| `mm_code_base_mismatch` | warning | `base_mismatch` | The mod code's canonical base differs from the base the segment declares |
| `mm_code_base_mismatch` | warning | `nonzero_on_inferred` | A position an earlier MM segment implied to be unmodified later receives a non-zero probability |

Two of the three `mm_semantics` checks are warnings: a read whose only semantic finding is `mm_code_unknown` or `mm_code_base_mismatch` still parses, its code is still counted, and `passed` stays true. `mm_position_unresolved` is the layer's one error, so a delta list that runs past the end of SEQ fails the run.

### `ml` — the ML probability array

| check_id | severity | anomaly_type | Fires when |
|----------|----------|--------------|------------|
| `ml_type_invalid` | error | `wrong_type` | ML is present but is not a `B:C` (byte array) aux field |
| `ml_too_short` | error | `too_short` | The ML array holds fewer values than MM declares |
| `ml_too_long` | error | `too_long` | The ML array holds more values than MM declares |

One segment declares `(number of mod codes) × (number of positions)` probabilities; the declared total is that summed over all segments.

## Check pipeline

Layers run in a fixed order for every record. A layer stops the ones below it only when its output is unusable -- not merely when it reports something:

| Layer | Runs when | On an unusable result |
|-------|-----------|-----------------------|
| `tag_set` | always | never stops a later layer |
| `mm_grammar` | the MM tag is present and is a string | no segments -> `mm_semantics` and `ml` do not run |
| `mm_semantics` | the MM string parsed | an incomplete site table -> `ml` does not run |
| `ml` | ML is present and the site table is complete | last layer |

What follows from this:

- A layer that never ran is reported exactly like one that ran and found nothing: `status: "ok"` with `issue_count: 0`. The report answers "what was caught", not "what was covered".
- A read can report findings from several layers at once -- the tag-set layer never blocks anything below it, and the semantic layer's two warnings are independent of the ML layer.

MM and ML pair as follows:

- **ML without MM** — `mm_missing` / `ml_without_mm`. The ML array itself is never inspected, so the `ml` layer stays `ok` however malformed that array is.
- **MM without ML** — `ml_missing` / `mm_without_ml`. The `ml` layer does not run, so this does not also report `ml_too_short` for the absent bytes.
- **Both present** — the ML byte count must equal the declared total.

## Report semantics

The same check runs automatically before every modbam-input command
(`convert`, `align`, `pileup`, `extract`, `summarize`, `evaluate`); see the
`Input validation` section of each command's page. That gate uses the
default 10000-record scan and never `--all`.

### Text mode

```
passed: True|False
reads_checked: <N>
mod_code_counts: a=50 m=800
scan: first 10000 records (use --all for the entire file)
<layer title>: ok | <E> error | <W> warning | <E> error, <W> warning
    <check_id>  <severity>  <issue_count>
      <anomaly_type>  count=<count>
        ERROR|WARNING   <read_id>  <tag>  <reason>
```

The four layer titles are `tag_set`, `MM grammar`, `MM semantics` and `ML`, printed in that order and always printed -- a clean file lists all four as `ok`. There is no blank line between the header lines and the first layer line.

- `mod_code_counts`: keys sorted lexicographically, printed as space-joined `code=count`; with no mod codes, only `mod_code_counts:` is printed, with no trailing token
- `scan:`: printed after `mod_code_counts` when `--all` is not set; omitted when `--all` is set
- check rows (`4` spaces): printed only for checks with `issue_count > 0`, in registry order within the layer
- anomaly-type rows (`6` spaces): sorted lexicographically
- sample rows (`8` spaces): the `ERROR` / `WARNING` prefix matches the check's `severity`; at most 5 rows per anomaly type, while `count` is the full count
- the separator inside a sample row is three spaces between the severity label and the read id, and two spaces elsewhere

### JSON mode (`--json`)

```json
{
    "passed": false,
    "reads_checked": 1,
    "mod_code_counts": {
        "m": 1
    },
    "layers": [
        {
            "id": "tag_set",
            "title": "tag_set",
            "status": "ok",
            "checks": {
                "mm_missing": {
                    "status": "ok",
                    "severity": "error",
                    "issue_count": 0
                },
                "ml_missing": {
                    "status": "ok",
                    "severity": "error",
                    "issue_count": 0
                },
                "tag_duplicate": {
                    "status": "ok",
                    "severity": "error",
                    "issue_count": 0
                },
                "tag_alias_conflict": {
                    "status": "ok",
                    "severity": "error",
                    "issue_count": 0
                }
            }
        },
        {
            "id": "mm_semantics",
            "title": "MM semantics",
            "status": "issues",
            "checks": {
                "mm_position_unresolved": {
                    "status": "issues",
                    "severity": "error",
                    "issue_count": 1,
                    "anomaly_types": {
                        "position_unresolved": {
                            "count": 1,
                            "samples": [
                                {
                                    "read_id": "r1",
                                    "tag": "MM",
                                    "reason": "invalid modification tag `MM`: delta list refers past end of sequence"
                                }
                            ]
                        }
                    }
                },
                "mm_code_unknown": {
                    "status": "ok",
                    "severity": "warning",
                    "issue_count": 0
                },
                "mm_code_base_mismatch": {
                    "status": "ok",
                    "severity": "warning",
                    "issue_count": 0
                }
            }
        }
    ]
}
```

The block above is abbreviated to two layers to show both states; a real run carries all four layers, and every check of every layer. `anomaly_types` is absent from a check with `issue_count == 0` and present otherwise. The tool prints the report with 4-space indentation, one key per line, in exactly the layout shown here.

**BREAKING:** `--json` no longer emits a flat `checks` array, and the top-level `errors` / `warnings` lists were removed in an earlier change. See [Migrating from the old report](#migrating-from-the-old-report).

### Field semantics

- **Nesting.** `layers` is an ordered array of the four layers, each with a stable machine-readable `id`, a display `title`, a `status`, and a `checks` mapping keyed by check id. Layer order is `tag_set`, `mm_grammar`, `mm_semantics`, `ml`; the `checks` mapping is in registry order, not lexicographic.
- **`status` is two-valued.** `status: "issues"` if and only if `issue_count > 0`. Both error and warning findings count towards `issue_count` and set `status` to `"issues"`; the level is expressed only by the sibling `severity` field. A layer has no `severity` of its own.
- **`severity` is a fixed property of the check**, taken from the tables under [Checks](#checks) rather than from the first finding. It is never decided per record.
- **The inventory is pre-filled.** Every layer and every check is present in every report with `status: "ok"`, `severity` set and `issue_count: 0` unless something fired. There is no field saying whether a check ran.
- **`passed`** is `false` if and only if some check with `severity: "error"` has `issue_count > 0`. Warnings never affect it, so `passed: true` coexists with a check reporting `status: "issues"`.
- **Bucketing.** Issues are aggregated by `(check_id, anomaly_type)`. `count` is the full number of issues of that subtype across the scanned reads; `samples` keeps at most **5** of them, first-come in scan order.
- **`mod_code_counts`** counts the mod codes declared by reads whose MM string parsed -- once per read, however many segments declare the code, and including codes whose probabilities are all zero. Unknown codes are counted too. Reads with no MM tag, reads whose MM is not a string, and reads whose MM failed to parse contribute nothing. Keys are lexicographic. The field does not affect `passed`.
- **Scan window.** Only the first 10000 records are checked unless `--all` is set. Records outside the window count towards nothing: not `reads_checked`, not `mod_code_counts`, not any check.

## Checks in detail

Each example starts from the same one-line SAM record and changes a single tag. Save it as `base.sam`:

```
@HD	VN:1.6	SO:coordinate
@SQ	SN:chr1	LN:1000
r1	0	chr1	1	60	16M	*	0	0	ACGTACGTACGTACGT	IIIIIIIIIIIIIIII	MM:Z:C+m,0,0;	ML:B:C,128,200
```

The read is 16 bp and `C` occurs at 0-based positions 1, 5, 9 and 13. `MM:Z:C+m,0,0;` declares two modified C's -- the first two, at positions 1 and 5 -- so MM requires exactly two ML values.

```bash
samtools view -b base.sam > base.bam
poremod validate --input base.bam
```

```
passed: True
reads_checked: 1
mod_code_counts: m=1
scan: first 10000 records (use --all for the entire file)
tag_set: ok
MM grammar: ok
MM semantics: ok
ML: ok
```

Every block below is quoted from a real run on this record. Only the header lines and the layer blocks that carry the example are reproduced; the `passed` / `reads_checked` / `mod_code_counts` / `scan` lines are as in the baseline unless a block includes them. A `mm_grammar` failure aborts before any modification code is recorded, so those examples show `mod_code_counts:` with no token.

### `mm_missing`

**`wrong_type`** — MM is not a string aux field. The check is decided on MM alone, so it fires whether or not ML is present; dropping the ML tag from the record below leaves the block byte-identical:

```
MM:i:5	ML:B:C,128,200
```

```
tag_set: 1 error
    mm_missing  error  1
      wrong_type  count=1
        ERROR   r1  MM  MM tag is not a string
```

Because no usable MM string exists, the `ml` layer never runs and `ML: ok` even though this read's ML array is a perfectly good two-byte one.

**`ml_without_mm`** — ML present, no MM tag. An orphaned ML is an error:

```
ML:B:C,128,200
```

```
tag_set: 1 error
    mm_missing  error  1
      ml_without_mm  count=1
        ERROR   r1  MM  ML tag present without an MM tag
```

### `ml_missing`

**`mm_without_ml`** — a usable MM string present, no ML tag:

```
MM:Z:C+m,0,0;
```

```
tag_set: 1 error
    ml_missing  error  1
      mm_without_ml  count=1
        ERROR   r1  ML  MM tag present without an ML tag
```

`mod_code_counts: m=1` — the MM string parsed, so its code is recorded even though the probabilities are missing.

### `tag_duplicate`

**`duplicate_mm`** — the read carries two MM tags. A SAM line may repeat a tag; `samtools` writes both through:

```
MM:Z:C+m,0,0;	MM:Z:C+m,0,0;	ML:B:C,128,200
```

```
tag_set: 1 error
    tag_duplicate  error  1
      duplicate_mm  count=1
        ERROR   r1  MM  read carries 2 MM tags
```

**`duplicate_ml`** — likewise for ML:

```
MM:Z:C+m,0,0;	ML:B:C,128,200	ML:B:C,128,200
```

```
tag_set: 1 error
    tag_duplicate  error  1
      duplicate_ml  count=1
        ERROR   r1  ML  read carries 2 ML tags
```

In both cases the first matching tag is the one read, so the layers below still run on a usable value and `mod_code_counts: m=1`.

### `tag_alias_conflict`

**`alias_mm`** — both `MM` and its lowercase alias `Mm`:

```
MM:Z:C+m,0,0;	Mm:Z:C+m,0,0;	ML:B:C,128,200
```

```
tag_set: 1 error
    tag_alias_conflict  error  1
      alias_mm  count=1
        ERROR   r1  MM  read carries both `MM` and `Mm`
```

**`alias_ml`** — both `ML` and `Ml`:

```
MM:Z:C+m,0,0;	ML:B:C,128,200	Ml:B:C,128,200
```

```
tag_set: 1 error
    tag_alias_conflict  error  1
      alias_ml  count=1
        ERROR   r1  ML  read carries both `ML` and `Ml`
```

A single alias on its own is legal: `MM:Z:C+m,0,0;` written as `Mm:Z:C+m,0,0;` passes silently.

### `mm_base_invalid`

**`invalid_base`** — the leading character must be one of `ACGTUN`:

```
MM:Z:X+m,0;	ML:B:C,128
```

```
MM grammar: 1 error
    mm_base_invalid  error  1
      invalid_base  count=1
        ERROR   r1  MM  MM segment `X+m,0`: base must be one of ACGTUN, got `X`
```

A segment shorter than two characters takes the same `check_id` and `anomaly_type`, with the reason `MM segment \`X\` is too short: expected a base and a strand`.

### `mm_strand_invalid`

**`invalid_strand`** — the second character must be `+` or `-`:

```
MM:Z:C*m,0;	ML:B:C,128
```

```
MM grammar: 1 error
    mm_strand_invalid  error  1
      invalid_strand  count=1
        ERROR   r1  MM  MM segment `C*m,0`: strand must be `+` or `-`, got `*`
```

### `mm_code_invalid`

**`invalid_code`** — three shapes share this anomaly type. A code that mixes letters with digits:

```
MM:Z:C+m0,0;	ML:B:C,128
```

```
MM grammar: 1 error
    mm_code_invalid  error  1
      invalid_code  count=1
        ERROR   r1  MM  MM segment `C+m0,0`: code `m0` must be all letters or all digits
```

An empty code, written `C+?,0;` or `C+,0;`:

```
MM:Z:C+?,0;	ML:B:C,128
```

```
MM grammar: 1 error
    mm_code_invalid  error  1
      invalid_code  count=1
        ERROR   r1  MM  MM segment `C+?,0`: modification code is empty
```

And a ChEBI code too large for `u32`, here `99999999999`, whose reason reads `MM segment \`C+99999999999,0\`: ChEBI code \`99999999999\` is out of range`.

### `mm_marker_invalid`

**`invalid_marker`** — `?` and `.` are the skip markers and may only be the last character of the code part:

```
MM:Z:C+?m,0;	ML:B:C,128
```

```
MM grammar: 1 error
    mm_marker_invalid  error  1
      invalid_marker  count=1
        ERROR   r1  MM  MM segment `C+?m,0`: `?`/`.` must be the last character of the code part
```

As the last character they are fine: `MM:Z:C+m?,0,0;` and `MM:Z:C+m.,0,0;` both pass.

### `mm_no_terminator`

**`missing_terminator`** — the MM string must end with `;`:

```
MM:Z:C+m,0,0	ML:B:C,128,200
```

```
MM grammar: 1 error
    mm_no_terminator  error  1
      missing_terminator  count=1
        ERROR   r1  MM  MM tag does not end with `;`
```

### `mm_empty_segment`

**`empty_segment`** — every segment carries exactly one terminator, so `;;` is an empty segment rather than a harmless extra:

```
MM:Z:C+m,0,0;;	ML:B:C,128,200
```

```
MM grammar: 1 error
    mm_empty_segment  error  1
      empty_segment  count=1
        ERROR   r1  MM  MM contains an empty segment
```

### `mm_delta_invalid`

**`invalid_delta`** — a delta list is `(,[0-9]+)*`, so every comma must be followed by digits. A leading comma is malformed rather than an absent delta:

```
MM:Z:C+m,,0;	ML:B:C,128
```

```
MM grammar: 1 error
    mm_delta_invalid  error  1
      invalid_delta  count=1
        ERROR   r1  MM  MM segment `C+m,,0`: delta `` is not a non-negative integer
```

A non-numeric delta hits the same anomaly type, with the reason `MM segment \`C+m,x\`: delta \`x\` is not a non-negative integer`. `C+m;` -- no comma at all -- legitimately declares zero deltas and passes the grammar layer.

### `mm_position_unresolved`

**`position_unresolved`** — the MM string parses, but a delta runs off the end of SEQ. The read holds four C's, and a first delta of 9 addresses the tenth:

```
MM:Z:C+m,9;	ML:B:C,128
```

```
MM semantics: 1 error
    mm_position_unresolved  error  1
      position_unresolved  count=1
        ERROR   r1  MM  invalid modification tag `MM`: delta list refers past end of sequence
```

This is the one semantic-layer error. It leaves the site table incomplete, so the `ml` layer does not run: note `mod_code_counts: m=1` alongside `ML: ok` despite the one-byte ML array.

### `mm_code_unknown`

**`unsupported`** — the mod code is not in the known-code table. This is a warning: the read still parses, and the code is still counted.

```
MM:Z:C+99999,0,0;	ML:B:C,128,200
```

```
passed: True
reads_checked: 1
mod_code_counts: 99999=1
scan: first 10000 records (use --all for the entire file)
tag_set: ok
MM grammar: ok
MM semantics: 1 warning
    mm_code_unknown  warning  1
      unsupported  count=1
        WARNING   r1  MM  unknown modification code '99999'
ML: ok
```

### `mm_code_base_mismatch`

**`base_mismatch`** — the code's canonical base differs from the base the segment declares. `m` is a C modification, so declaring it on an `A` segment warns:

```
MM:Z:A+m,0;	ML:B:C,128
```

```
MM semantics: 1 warning
    mm_code_base_mismatch  warning  1
      base_mismatch  count=1
        WARNING   r1  MM  modification code 'm' is a C modification but the segment declares base A
```

A segment declared on `N` is a wildcard and never mismatches. ChEBI codes and the four "any" codes have no fixed base, so they are never reported here.

**`nonzero_on_inferred`** — an MM segment that carries no explicit skip marker implicitly declares every canonical base it does not list as unmodified. A later segment of the same base and strand must not then assign a non-zero probability to one of those positions:

```
MM:Z:C+m,0,0;C+m,2;	ML:B:C,0,0,128
```

```
passed: True
reads_checked: 1
mod_code_counts: m=1
scan: first 10000 records (use --all for the entire file)
tag_set: ok
MM grammar: ok
MM semantics: 1 warning
    mm_code_base_mismatch  warning  1
      nonzero_on_inferred  count=1
        WARNING   r1  MM  implicit canonical position 9 (+) carries a non-zero modification probability
ML: ok
```

The first segment covers positions 1 and 5, and its implicit fill marks positions 9 and 13 as unmodified. The second segment then covers position 9 with a probability of 128 rather than 0. Warning only: `passed` stays `True` and the exit code is 0. Setting the third ML value to 0 removes the warning.

### `ml_type_invalid`

**`wrong_type`** — ML must be a byte array (`B:C`). Any other aux type, here a string:

```
MM:Z:C+m,0,0;	ML:Z:abc
```

```
ML: 1 error
    ml_type_invalid  error  1
      wrong_type  count=1
        ERROR   r1  ML  ML tag is not a u8 array
```

### `ml_too_short`

**`too_short`** — the ML array holds fewer bytes than MM declares. Truncate it to one byte:

```
MM:Z:C+m,0,0;	ML:B:C,128
```

```
ML: 1 error
    ml_too_short  error  1
      too_short  count=1
        ERROR   r1  ML  ML array too short: segment 1 (`C+m,0,0`) needs 2 probabilities, only 1 remain
```

The shortfall is attributed to the segment where the cumulative total first exceeds the array, and the reason names it by 1-based index. A segment shorter than 40 characters is spelled out in full; a longer one is cut at a delta boundary and closed with `...`, so a segment listing thousands of calls cannot flood the report.

### `ml_too_long`

**`too_long`** — extend the same array to three bytes:

```
MM:Z:C+m,0,0;	ML:B:C,128,200,255
```

```
ML: 1 error
    ml_too_long  error  1
      too_long  count=1
        ERROR   r1  ML  ML array too long: MM declares 2 probabilities, ML has 3
```

### Several layers at once

A read can be reported by more than one layer. The tag-set layer never blocks anything below it, so a record can carry a tag-pairing error and a grammar error at once -- here an omitted ML tag alongside an unparseable MM string:

```
MM:Z:C+m0,0;
```

```
tag_set: 1 error
    ml_missing  error  1
      mm_without_ml  count=1
        ERROR   r1  ML  MM tag present without an ML tag
MM grammar: 1 error
    mm_code_invalid  error  1
      invalid_code  count=1
        ERROR   r1  MM  MM segment `C+m0,0`: code `m0` must be all letters or all digits
```

A semantic-layer warning does not stop the ML layer, so combining an unknown code with a short array hits a warning and an error:

```
MM:Z:C+99999,0,0;	ML:B:C,128
```

```
MM semantics: 1 warning
    mm_code_unknown  warning  1
      unsupported  count=1
        WARNING   r1  MM  unknown modification code '99999'
ML: 1 error
    ml_too_short  error  1
      too_short  count=1
        ERROR   r1  ML  ML array too short: segment 1 (`C+99999,0,0`) needs 2 probabilities, only 1 remain
```

A layer header aggregates both levels when they appear together, as `MM semantics: 1 error, 1 warning`.

## Migrating from the old report

`poremod validate --json` no longer emits a flat `checks` array. Reads that used `report["checks"]` must switch to walking `report["layers"]`:

```python
for layer in report["layers"]:
    for check_id, check in layer["checks"].items():
        if check["issue_count"]:
            ...
```

`check_id -> severity` is now a fixed mapping (see the tables under [Checks](#checks)) rather than a value taken from the first finding. `checks` is a mapping keyed by check id, not a list; the flat per-check `severity` / `issue_count` / `anomaly_types` fields survive inside each entry — the check id is no longer a field but the key of the enclosing `checks` mapping — while `anomaly_types` is itself a mapping keyed by anomaly type rather than a list of `{"type": ...}` objects. The top-level `errors` / `warnings` lists were removed earlier and are still absent.

Some ids were renamed outright, and a script that indexes a check by its old id has to be updated even though nothing else about it changed:

| Old id (anomaly type) | New id (anomaly type) |
|-----------------------|-----------------------|
| `ml_length` (`too_short` / `too_long`) | `ml_too_short` (`too_short`) / `ml_too_long` (`too_long`) |
| `unknown_mod_code` (`unsupported`) | `mm_code_unknown` (`unsupported`) |
| `implicit_canonical` (`nonzero_on_inferred`) | `mm_code_base_mismatch` (`nonzero_on_inferred`) |
| `mm_syntax` (`parse_error`) | `mm_grammar` splits this coarse id across its own: `mm_base_invalid`, `mm_strand_invalid`, `mm_code_invalid`, `mm_marker_invalid`, `mm_no_terminator`, `mm_empty_segment`, `mm_delta_invalid`; the "MM is present but not a string" case is now `mm_missing` (`wrong_type`) |
| `mm_syntax` (`position_error`) | `mm_position_unresolved` (`position_unresolved`) |
| `ml_type` (`wrong_type`) | `ml_type_invalid` (`wrong_type`) |

Two items were retired rather than renamed. **`mn_tag` is gone: MN is no longer checked at all**, by no layer, so a pipeline that relied on the old length check silently loses it. The `htslib_error` anomaly type of `mm_syntax` / `ml_type` (and of `mn_tag`) is gone too — an htslib failure is never reported as an anomaly in the layered report. One further old behaviour moved rather than disappeared: `ml_length` (`too_short`) also fired for a read whose ML tag was absent while MM declared probabilities, and that read is now the tag-set layer's `ml_missing` (`mm_without_ml`), with no `ml_too_short` emitted alongside it.

## Exit codes

| Exit code | Scenario |
|-----------|----------|
| 0 | No errors (warnings only / all passed) |
| 1 | At least one error-level check |
| 2 | Input missing / cannot be opened |

## Examples

```bash
poremod validate --input in.bam
poremod validate --input in.bam --all
poremod validate --input in.bam --json > report.json
```
