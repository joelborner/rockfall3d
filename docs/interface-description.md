# Executable Interface Description, Version 1

## Command and Coverage

```text
rockfall3d.exe --describe-interface
```

The command returns one UTF-8 JSON document on stdout, no byte-order mark, and
exit code 0 on success. Diagnostics belong on stderr. It requires no scenario
or input data, starts no simulation, and creates no files. Normal runtime DLL
dependencies still have to be installed beside the executable or discoverable.

The source is `interface/description.json.in`. CMake embeds it in the binary;
users do not need a JSON sidecar. The generated copy in the build directory is
for inspection only. Capability selection must describe actual compiled support,
never infer it from branch names or filenames of selected executables.

This pilot manages only **optional rock models**, domain `rock.models`. The
current experimental branch advertises one Fragmentation section, disabled by
default. Builds without fragmentation return an empty `sections` list. This does not disable existing density, scaling, release,
terrain, numerical or output controls: those are outside the managed domain.

## Document

| Property | Meaning |
| --- | --- |
| `interface_schema_version` | Integer `1`; the structure/rules version, not the simulator version. |
| `engine.id` | Stable simulator ID, `rockfall3d`. |
| `engine.version` | Software version from the CMake project. |
| `engine.build_id` | Git revision, suffixed `-dirty` for modified/untracked source; `source-archive` when unavailable. |
| `engine.source_revision` | Git revision string, or `null` when unavailable. |
| `engine.source_dirty` | Boolean, or `null` when unavailable. |
| `engine.description_id` | SHA-256 of the description source, prefixed `sha256:`. Fixtures use explicit fixture IDs. |
| `managed_domains` | Domains for which this description gives a complete supported-section catalogue. |
| `sections` | Supported sections in those domains. |

Within a managed domain, an absent section is unsupported. Outside those
domains, absence means not described yet. Two dirty builds can share a build ID;
it is provenance, not a unique executable fingerprint. Cache by executable
identity (path plus replacement detection), not build ID alone. Source archives
also cannot establish exact provenance. The description hash identifies the
catalogue source, independently of Git revision.

Provenance refreshes at build time. An unchanged build does not gain new identity
merely because it is queried. The description is never loaded from the source
tree at runtime.

## Sections and Parameters

Each section has:

- `id`: stable machine ID, unique across the document.
- `domain`: one of `managed_domains`.
- `target`: a registered semantic destination; the pilot uses `rock`.
- `title` and `description`: readable text/help.
- `order`: integer ordering within its destination, then ID as a tie breaker.
- `optional`: whether the user can disable the section.
- `enabled_by_default`: initial enabled state for a new scenario.
- `parameters`: ordered field definitions.
- `rules`: cross-parameter validation rules, all of which must pass when enabled.
- `xml`: section mapping, described below.

The plugin owns tab names, placement relative to specialised widgets, and layout.
`target` adds a section to an existing registered destination; it does not create
arbitrary application tabs. A future terrain-card container can reuse the form
renderer for each property, but scopes, inheritance and spatial processing are
not part of this pilot.

Each parameter has:

| Property | Meaning |
| --- | --- |
| `id` | Stable ID, unique within the section; store values by section ID and parameter ID. |
| `label`, `description` | Readable label and help. |
| `type` | The initial v1 vocabulary supports `number` (finite real values, excluding booleans). |
| `unit` | Display unit; `1` denotes dimensionless and may be displayed as `[-]`. |
| `initial_value` | Starting GUI value for a new section, never a substitute for missing required imported XML. |
| `required` | Whether the XML attribute must be present whenever the element is present. |
| `parser_default` | Effective engine value when an optional XML attribute is omitted; absent for required attributes. |
| `minimum`, `maximum` | Optional inclusive bounds. |
| `exclusiveMinimum`, `exclusiveMaximum` | Optional exclusive bounds. |
| `xml.attribute` | Accepted XML attribute name, independent of the readable label and parameter ID. |

All declared bounds apply. Avoid arbitrary widget precision or range limits that
silently clamp valid values. Editing/display precision must preserve values on
save and reopen. Rules are evaluated only for enabled sections. For an existing
XML element, a missing required attribute is an error; optional attributes use
`parser_default`. The writer can explicitly emit all resolved parameter values.

The numeric names follow JSON Schema conventions; this document is a small
application contract, not a JSON Schema validator or an executable UI language.
Additional field types and rule operators require an agreed contract extension
and plugin support; never guess their behaviour.

## Structured Validation

The pilot defines exactly one cross-parameter operator:

```json
{
  "operator": "greater_than",
  "left": {"parameter": "saturation_specific_energy"},
  "right": {"parameter": "threshold_specific_energy"},
  "message": "Saturation specific energy must exceed threshold specific energy."
}
```

Operands reference parameters in the same section, after optional defaults are
resolved. The rule means strict numeric `left > right`. Use the message when it
fails. There are no expression strings, scripts, templates, or arbitrary code to
evaluate. Unknown operators must not be silently skipped.

## XML Mapping

```json
{
  "parent_path": [],
  "element": "fragmentation",
  "emit_when": "enabled",
  "cardinality": "zero_or_one"
}
```

`parent_path` lists enclosing element names from the document level. The empty
path means a document-level element. In the existing simulator format,
`<rockfall3d/>` is a marker followed by sibling configuration elements, not an
enclosing root. Do not put fragmentation inside that self-closing marker. A
reader using an XML library may need a synthetic root to read the legacy
fragment format; it must not change the emitted simulation syntax.

`emit_when: enabled` omits the element completely when disabled. When enabled,
emit one self-closing element and its mapped numeric attributes. Attribute order
and whitespace are not significant. `zero_or_one` forbids duplicate model
elements. During import, an existing element enables its section. Use an XML
library to escape and serialize attributes rather than evaluating templates.

## Executable Switching and Legacy Builds

Query the selected executable asynchronously, with a timeout, and keep stdout
and stderr separate. Validate the document before using it. An unknown schema,
field type, rule or target requires a clear compatibility message; do not create
a partially understood runnable section.

Only sections advertised by the selected executable appear, are validated, and
produce model XML. Keep other stored values internally for switching back; they
must not resurrect unsupported forms, affect current validation, or block a run.
The optional section's checkbox controls activation only after support has been
established from the description. Old feature-specific preference flags do not
establish capabilities. Values outside managed domains remain the responsibility
of the existing plugin controls. Save the description snapshot and engine
identity with scenario settings for traceability.

The plugin requires a successful description query before previewing/staging or
running a simulation. An older executable without this command has unknown
capabilities and cannot use a saved or bundled catalogue as a runtime fallback.
A failed query can also mean a missing runtime DLL, timeout or invalid executable.
Bundled definitions may interpret historical XML, but must not decide which
controls the selected executable provides.

## Fragmentation Pilot Fixture and Integration

Complete fixtures live only in the companion `rockfall3d-qgis` repository:

- `tests/fixtures/interface/without-fragmentation.json`: a build without optional rock models.
- `tests/fixtures/interface/with-fragmentation.json`: the supported experimental section.
- `tests/fixtures/interface/fragmentation.xml`: corresponding example element.

They test the plugin's actual description validator, generated forms, executable
switching and XML mapping. They are not input to C++ tests or to production
capability discovery. The fragmentation fixture describes these existing
experimental-branch inputs:

| ID | XML attribute | Initial value | Parser default | Constraint |
| --- | --- | --- | --- | --- |
| `threshold_specific_energy` | `threshold_specific_energy` | 0.2 J/kg | Required | >= 0 |
| `saturation_specific_energy` | `half_specific_energy` | 100 J/kg | Required | > threshold |
| `exponent` | `exponent` | 2 | 2 | > 0 |
| `max_removed_fraction` | `max_removed_fraction` | 0.5 | 0.5 | > 0 and <= 0.5 |

`half_specific_energy` is retained for compatibility; the model reaches maximum
removal there. The GUI label deliberately says saturation. Do not expose the
shape directory; a newly generated element omits `shape_output_dir` and uses
the engine default. On import, preserve any explicitly supplied hidden directory
attribute when saving an otherwise unchanged scenario. Preserve the engine's
automatic mass/shape-index output behaviour.

The fixture is **not a capability switch**. This branch registers the section in
its actual `interface/description.json.in` and compiles the fragmentation model.
The built executable therefore supplies the section directly; no plugin override
or external JSON file is needed.

Metadata and parser definitions remain separate. `fragmentation_interface_tests`
feeds the advertised XML names, initial values and defaults into tests of the real
parser. It also checks required attributes, numeric limits, omission, the default
shape directory and automatic mass/shape-index outputs. These are engine tests,
not hypothetical plugin fixtures. A future refactor can share definitions more
directly without rewriting the full parser for this pilot.

## Verification

Build normally, then run:

```powershell
ctest --test-dir build/windows-release -C Release --output-on-failure
```

`interface_description_tests` uses CMake to check the real executable's JSON,
identity, clean stdout/stderr, operation from an empty working directory, absence
of generated files, argument errors and agreement with the generated catalogue.
It requires no Python, plugin imports, or model fixtures. It does not assume
that fragmentation, or any other optional section, exists in a particular build.

The plugin project owns the fixture-based tests for numeric bounds, nonfinite
values, required/default semantics, cross-parameter rules, optional-section
omission, generated forms, executable switching and XML round trips. These tests
exercise production plugin code rather than a separate reference implementation.

The end-to-end plugin acceptance check is: switch builds, show/hide generated
fragmentation, retain stored values, save/reopen XML, and run accepted input with
this branch's executable. Existing simulation invocation,
numerical integration, solvers and output formats are unchanged by this command.
