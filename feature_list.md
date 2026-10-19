# GWCS.jl feature checklist

This document tracks useful Python GWCS capabilities that GWCS.jl does not yet support. It also includes a small number of Julia-oriented interoperability features that do not have a direct Python equivalent but are needed to make the package broadly useful.

Priorities:

- **P0**: required to read common real-world GWCS metadata or preserve correctness
- **P1**: broadly useful core functionality
- **P2**: specialized, optional, or suitable for a package extension

Python reference surfaces include `gwcs.coordinate_frames`, `gwcs.wcs`, `gwcs.wcstools`, `gwcs.selector`, `gwcs.geometry`, `gwcs.spectroscopy`, `gwcs.fitswcs`, and the converters in `gwcs.converters`.

## Transformation protocol and composition

- [x] **P0** Add a public transform-registration API so external packages can register metadata tags, constructors, inverses, dimensions, and compilation rules without modifying GWCS.jl.
- [x] **P1** Add a callable transform wrapper for user functions with declared input/output dimensions and optional inverse and Jacobian functions.
- [x] **P1** Add constant-output and fixed-input transforms, including the ASDF `constant` and `fix_inputs` semantics.
- [ ] **P1** Extend axis remapping to support dropping inputs, inserting constants, and named input/output axes.
- [x] **P1** Support elementwise arithmetic composition of compatible transforms: add, subtract, multiply, divide, and power.
- [x] **P1** Preserve transform names and input/output axis names through metadata parsing, composition, inversion, and compilation.
- [x] **P1** Add per-transform dependency or separability metadata and compute an axis-correlation matrix for a compiled pipeline.
- [x] **P2** Add a stable canonical tagged transform-tree inspection API for debugging and serialization.
- [ ] **P2** Add simplification passes beyond affine fusion, such as identity removal across nested transforms, constant folding, remap fusion, and polynomial/linear fusion where numerically safe.

## Celestial projections and rotations

- [x] **P0** Parse ASDF `gnomonic` projection nodes directly into the existing TAN transforms.
- [x] **P0** Parse ASDF `rotate3d` nodes for native-to-celestial, celestial-to-native, and general Euler rotations.
- [x] **P0** Parse spherical rotation sequences in addition to the currently supported Cartesian rotation sequences.
- [x] **P1** Add the remaining FITS Paper II zenithal projections: AZP, SZP, STG, SIN, ARC, ZEA, and AIR.
- [x] **P1** Add cylindrical projections: CYP, CEA, CAR, and MER.
- [x] **P1** Add pseudocylindrical projections: SFL, PAR, MOL, and AIT.
- [x] **P1** Add conic projections: COP, COE, COD, COO, and BON.
- [x] **P1** Add polyconic and spherical-cube projections: PCO, TSC, CSC, and QSC.
- [x] **P1** Add HEALPix-family projections: HPX and XPH.
- [x] **P1** Support zenithal projection parameters and explicit inverses with the same domain conventions as Astropy and wcslib.
- [x] **P1** Support parameters and explicit inverses for the remaining projection families with the same domain conventions as FITSWCS.jl and wcslib.
- [ ] **P1** Add a general 2D rotation transform and named Euler-angle rotation transform.
- [x] **P1** Add automatic `LONPOLE` computation and general native-pole solving for supported celestial projections.
- [ ] **P1** Support celestial reference frames beyond ICRS, including Galactic, FK4, FK5, ecliptic, and supergalactic frames.
- [ ] **P2** Add optional conversion between Julia celestial-coordinate objects and numeric longitude/latitude pipeline values.

## Polynomial, interpolation, and mathematical transforms

- [x] **P0** Parse ASDF `affine` nodes, including explicit translation and inverse nodes, using the existing `AffineTransform`.
- [x] **P0** Parse ASDF `identity` nodes.
- [x] **P0** Parse legacy and current ASDF `constant` nodes.
- [ ] **P1** Add Chebyshev polynomials in one and two dimensions with domain/window mapping.
- [ ] **P1** Add Legendre polynomials in one and two dimensions with domain/window mapping.
- [ ] **P2** Add Hermite polynomials in one and two dimensions.
- [ ] **P1** Add one- and two-dimensional tabular interpolation transforms with configurable interpolation method, bounds behavior, and fill value.
- [ ] **P2** Generalize tabular interpolation to arbitrary dimensions where a concrete use case requires it.
- [ ] **P2** Add one-dimensional spline transforms and their ASDF representation.
- [ ] **P1** Add serializable elementary mathematical transforms used in metadata-defined pipelines: trigonometric, inverse trigonometric, logarithmic, exponential, square root, absolute value, reciprocal, and degree/radian conversion.
- [x] **P1** Add a named SIP distortion constructor and inverse representation, even if it compiles internally to shifts and polynomials.
- [x] **P2** Add native TPV and Paper IV TPD distortion transforms.

## Geometry, spectroscopy, and optical transforms

- [ ] **P1** Add `ToDirectionCosines` and `FromDirectionCosines`.
- [ ] **P2** Add the grating-equation transforms `WavelengthFromGratingEquation` and `AnglesFromGratingEquation3D`.
- [ ] **P2** Add the `Snell3D` refraction transform.
- [ ] **P2** Add the `SellmeierGlass` dispersion transform.
- [ ] **P2** Add the temperature- and pressure-dependent `SellmeierZemax` transform.
- [ ] **P2** Place specialized spectrograph and optical transforms in an optional extension or companion package unless they become common dependencies of multiple observatory pipelines.

## Discontinuous and region-selected transforms

- [x] **P1** Add a region-selector transform that dispatches to different compiled transforms based on a label.
- [x] **P1** Add array-backed label maps for detector masks and IFU slices.
- [x] **P1** Add dictionary-backed label maps with configurable numeric matching tolerance.
- [x] **P1** Add range-backed label maps with overlap validation.
- [x] **P1** Add function- or transform-backed label maps.
- [x] **P1** Support a configurable output value when no region transform is defined.
- [x] **P1** Construct an inverse region selector when the label mapper and every selected transform provide compatible inverses.
- [x] **P1** Add selector-dependent compound bounding boxes and ignored bounding-box inputs.
- [ ] **P2** Add polygon-to-label-mask utilities using a robust geometry implementation rather than copying Python GWCS's scan-line implementation.
- [ ] **P2** Add helpers for selecting or extracting the transform associated with one IFU slice, slit, spectral order, or other region label.

## Coordinate frames, units, and high-level coordinates

- [x] **P0** Make `axes_order` operational during evaluation and inverse evaluation rather than storing it only as metadata.
- [x] **P1** Add an `AxisType` representation for spatial, spectral, temporal, Stokes, pixel, undefined, and custom axes.
- [x] **P1** Add `SpectralFrame` with wavelength, frequency, energy, and velocity physical-type inference.
- [x] **P1** Add `TemporalFrame` with a reference epoch, time scale, location metadata, and offset conversion.
- [x] **P1** Add `StokesFrame`.
- [x] **P1** Add `CompositeFrame` for mixed outputs such as celestial plus spectral or celestial plus spectral plus time.
- [x] **P1** Validate that composite-frame axis orders are unique and complete.
- [x] **P1** Add public accessors for pixel/world axis names, physical types, units, and axis types.
- [ ] **P1** Integrate frame-boundary unit validation and conversion through an optional Unitful.jl extension.
- [ ] **P2** Integrate celestial outputs with an established Julia coordinate package through an optional extension.
- [ ] **P2** Integrate temporal outputs with an established Julia time/astronomy package through an optional extension.
- [ ] **P2** Define an idiomatic Julia protocol for converting raw numeric axes to and from richer world-coordinate objects.
- [ ] **P2** Support non-coordinate inputs and outputs, such as spectral-order selectors, without pretending that they are physical coordinate axes.
- [ ] **P2** Validate physical-type strings against an IVOA UCD vocabulary when validation data is available.

## Evaluation APIs and WCS introspection

- [x] **P1** Add `get_frame(wcs, name)` as a public frame lookup.
- [ ] **P1** Add `pixel_to_world_values` and `world_to_pixel_values` aliases for explicitly numeric results.
- [ ] **P1** Add array-index APIs that reverse Cartesian/array axis order and round world-to-index results consistently.
- [ ] **P1** Add broadcasted coordinate-array evaluation where each axis is supplied as a scalar or array with compatible shapes.
- [ ] **P1** Add allocation-conscious in-place batch APIs for large coordinate arrays.
- [ ] **P1** Add `pixel_bounds` and an `array_shape` view of `pixel_shape`.
- [x] **P1** Add `axis_correlation_matrix(wcs)` based on transform dependency information.
- [x] **P1** Add `in_image(wcs, world...)` to test whether world coordinates invert to finite pixels inside the valid domain.
- [x] **P1** Add `in_domain(wcs, pixel...)` and matrix masks for finite coordinates inside the input domain.
- [x] **P1** Add `fill_invalid=true` so `pixel_to_world` and `world_to_pixel` return `fill_value` for coordinates whose transforms raise `DomainError`.
- [ ] **P1** Make invalid-coordinate filling fast: catching a thrown `DomainError` costs tens of microseconds per point (about 79 us per off-sky pixel for an AIT all-sky grid), which is impractical for all-sky maps and reprojection. Likely requires projections to report validity without throwing, selected when the WCS is constructed.
- [ ] **P1** Add explicit output-footprint masks after footprint and celestial-boundary semantics are available.
- [ ] **P1** Add `footprint(wcs)` for spatial, spectral, temporal, custom, or all output axes.
- [ ] **P1** Handle longitude wrapping and pole crossings when computing and testing celestial footprints.
- [x] **P2** Add concise text and MIME displays showing frames, transforms, dimensions, inverses, and bounding boxes.
- [ ] **P2** Add transform and pipeline visualization suitable for documentation and debugging.

## Numerical inversion

- [ ] **P1** Add a batched numerical-inverse kernel rather than invoking the scalar solver independently for every matrix column.
- [ ] **P1** Track convergence, slow convergence, divergence, achieved accuracy, and iteration count per point.
- [ ] **P1** Extend `NoConvergence` with per-point divergent and slow-convergence indices for batch calls.
- [ ] **P1** Add adaptive iteration that stops updating points that have already converged.
- [ ] **P1** Add explicit divergence detection and preserve the best pre-divergence solution.
- [ ] **P1** Fit and cache an approximate inverse transform for expensive distortion models when a bounding box is available.
- [ ] **P1** Use the approximate inverse as an initial guess before falling back to the bounding-box center.
- [ ] **P1** Allow a caller-supplied array of initial guesses for batch inversion.
- [ ] **P2** Add a solver abstraction so alternative Julia root-finding methods can be selected without changing the WCS API.
- [ ] **P2** Add optional Jacobian caching or analytical-Jacobian hooks for transforms that can provide them.
- [x] **P1** Reject world coordinates that the analytically invertible tail of a numerically inverted pipeline cannot reach (for example the TAN antipode of a TAN-SIP WCS) with `DomainError` before iterating.

## Pipeline construction and mutation

- [x] **P1** Add a checked pipeline builder that can be edited and then compiled into an immutable `WCS`.
- [x] **P1** Add `set_transform` for replacing the transform between adjacent frames.
- [x] **P1** Add `insert_transform` before or after an existing frame.
- [x] **P1** Add `insert_frame` for splitting an existing pipeline edge with a new intermediate frame.
- [x] **P1** Add removal of an intermediate frame while composing its adjacent transforms.
- [x] **P1** Add frame renaming with uniqueness checks.
- [x] **P1** Add append and prepend operations for frames and transforms.
- [x] **P1** Add persistent, non-mutating variants that return a newly compiled `WCS`.
- [x] **P1** Revalidate dimensions, frame uniqueness, terminal-step rules, bounds, and inverse availability after every edit.
- [x] **P1** Make recompilation atomic so a failed edit never leaves a partially valid WCS.
- [x] **P1** Add `fix_inputs` to bind selected transform inputs to constants and return a lower-dimensional WCS.
- [ ] **P2** Preserve unchanged compiled pipeline segments across edits where doing so measurably reduces construction cost.
- [ ] **P2** Add a pipeline diff or summary facility for auditing metadata-driven edits.

## Bounding boxes, domains, and sampling

- [x] **P1** Add compound bounding boxes keyed by selector values.
- [x] **P1** Support ignored inputs in bounding-box definitions.
- [x] **P1** Add public setters or reconstruction helpers for changing a WCS bounding box and pixel shape.
- [x] **P1** Add `grid_from_bounding_box` with scalar or per-axis step size and pixel-center handling.
- [ ] **P1** Support unit-bearing bounding-box intervals through the optional units extension.
- [ ] **P2** Add general domain objects beyond rectangular boxes, including polygonal and mask-backed domains.
- [ ] **P2** Add tools to transform or approximate a domain between intermediate frames.

## Programmatic WCS construction and fitting

- [ ] **P1** Add `wcs_from_fiducial` for constructing celestial, spectral, and composite WCS pipelines from a reference coordinate.
- [ ] **P1** Add a FITS-imaging-style constructor accepting `crpix`, `crval`, `cdelt`, `pc` or `cd`, projection, and pole parameters.
- [ ] **P1** Add `wcs_from_points` for fitting a celestial WCS from matched detector and sky coordinates.
- [ ] **P1** Support monomial, Chebyshev, and Legendre distortion fits in `wcs_from_points`.
- [ ] **P1** Fit and attach an approximate inverse during point-based WCS construction.
- [ ] **P1** Report fit residuals, conditioning, selected degree, and achieved forward/inverse error.
- [ ] **P2** Add robust fitting options for outlier rejection and weighted calibration points.
- [ ] **P2** Add constructors for common mixed celestial/spectral instrument geometries.

## ASDF and tagged metadata support

- [x] **P0** Parse ASDF tags for transforms already implemented programmatically: identity, affine, TAN/gnomonic, native/celestial rotation, and explicit 3D Euler rotation.
- [x] **P0** Parse all GWCS frame tags: generic frame, composite frame, spectral frame, temporal frame, and Stokes frame.
- [x] **P0** Parse region-selector and label-mapper GWCS tags.
- [ ] **P0** Parse direction-cosine, FITS-imaging, and spectroscopy GWCS tags as their transform implementations become available.
- [x] **P0** Parse ASDF fixed-input compound-transform tags.
- [x] **P0** Parse ASDF compound-transform tags for arithmetic operations.
- [x] **P0** Parse all implemented zenithal projection tags.
- [ ] **P0** Parse orthogonal-polynomial, tabular, spline, mathematical-function, and remaining projection tags as their implementations become available.
- [x] **P0** Support compound bounding-box metadata.
- [x] **P0** Validate tag properties and dimensions against supported schema versions rather than relying only on constructor errors.
- [x] **P0** Add explicit schema-version dispatch and migrations for known incompatible metadata versions.
- [ ] **P0** Preserve unknown optional metadata fields so a read/write cycle does not silently discard them.
- [x] **P1** Add canonical WCS, frame, transform, inverse, bounding-box, and pixel-shape serialization.
- [x] **P1** Add writing of canonical inline metadata through ASDF.jl.
- [ ] **P1** Add lazy, external, compressed, and explicit ASDF array-block writing policies.
- [x] **P1** Add canonical round-trip tests for the supported Python-produced YAML fixtures and ASDF.jl files.
- [ ] **P1** Add round-trip tests against Python-produced ASDF files for every supported tag.
- [ ] **P1** Add an extension API for observatory- or instrument-specific ASDF schemas.
- [ ] **P2** Add utilities to inspect an unsupported metadata tree and report all missing tags before parsing begins.

## FITS interoperability

- [x] **P0** Add native `parse_fits_wcs` and `find_all_fits_wcs` interfaces that construct ordinary compiled GWCS pipelines without a second execution hierarchy.
- [x] **P0** Parse FITS image-header axis counts, CRPIX/CRVAL/CDELT, PC/CD/CROTA on the celestial longitude/latitude block, WCSAXES-only primary and alternate keys, units, names, and partial image shapes.
- [x] **P0** Cover linear, celestial, linear spectral, linear time, Stokes, and mixed/interleaved axis WCSes.
- [x] **P1** Add native SIP, TPV, and Paper IV TPD transforms, including SCAMP TAN-plus-PV recognition and numerical inversion when explicit inverse coefficients are absent.
- [x] **P1** Add Paper IV D2IM/CPDIS lookup distortions, including NAXES/AXIS maps and one-dimensional lookup images, with construction-time external image-array resolution.
- [x] **P1** Add statically dimensioned multilinear `-TAB` transforms with direct separable monotonic 1D inverses and checked numerical inversion for coupled tables.
- [x] **P1** Add Paper III nonlinear spectral algorithms and grism models.
- [ ] **P1** Add `to_fits_sip` to fit a two-dimensional celestial GWCS with configurable forward and inverse degree/error limits.
- [ ] **P1** Add `to_fits_tab` for non-celestial and multidimensional WCS approximations.
- [ ] **P1** Add a combined `to_fits` path that chooses linear, SIP, and `-TAB` representations by separable axis group.
- [ ] **P1** Return approximation diagnostics and verify the generated FITS representation over the requested bounding box.
- [x] **P1** Preserve FITS 1-based pixel-origin conventions explicitly in the parsed affine stage.
- [x] **P2** Add optional FITSIO.jl and FITSFiles.jl adapters for detached headers, linked HDUs, and construction-time lookup-image resolution without making either package a core dependency.
- [x] **P2** Parse direct Paper IV SIP distortion records (`CPDISja='SIP'` with repeated `DPja` fields such as `SIP.FWD.p_q` and `SIP.REV.p_q`), normalizing them to the existing SIP forward and approximate-reverse representation.
- [x] **P2** Use `SIP.REV` and `TPD.REV` as inverse initializers, then refine against the forward distortion.
- [x] **P2** Import the remaining legacy polynomial distortions: DSS and WAT-encoded TNX/ZPX.
- [x] **P2** Support generic FITS logarithmic axes, including `TIME-LOG`, with a dedicated analytical `FITSLogTransform`.
- [x] **P2** Recognize WCSLIB-compatible generic celestial axis identifiers such as `XLON`/`XLAT`, including reordered axes and strict matching of `xLON`/`xLAT` or `xxLN`/`xxLT` pairs, while retaining richer metadata for known coordinate systems.
- [ ] **P3** Parse binary-table image-array and pixel-list WCS keyword families with explicit table-column selection; keep this separate from image-header parsing. Requires a dedicated table-column selection API and keyword grammar rather than extending image-header inference.
- [ ] **P3** Preserve complete FITS time-coordinate metadata, including `DATEREF`, `TIMEOFFS`, observation intervals, uncertainty/resolution fields, additional time scales, and geodetic observer locations. Requires designed time/location frame interfaces before this metadata can be exposed usefully.
- [x] **P2** Apply standard `RADESYS` and `EQUINOX` defaults for equatorial and ecliptic frames, including FK4/FK5 defaults and clearing inapplicable metadata for ICRS, GAPPT, Galactic, and non-celestial axes.
- [ ] **P1** Parse the FITS unit grammar rather than relying on enumerated spellings. Requires the planned optional unit integration or a deliberately scoped FITS-unit parser.
- [ ] **P3** Support the standard `CUBEFACE` axis representation for TSC, CSC, and QSC cube projections. Requires a multidimensional cube-face transform and explicit face-validity semantics.
- [ ] **P3** Preserve coordinate-error and phase metadata including `CRDERi`, `CSYERi`, `CZPHSi`, and `CPERIi`. Requires frame metadata interfaces for coordinate uncertainty and phase conventions.
- [ ] **P3** Support Paper IV's arbitrary-power general `Polynomial` distortion when TPD normalization is not exact.
- [ ] **P3** Support sequent `CQDIS=LOOKUP` distortions and construction-time resolution of their lookup arrays.
- [x] **P1** Add an explicit FITS compatibility/repair policy for legacy aliases and conventions, including `NCP`, `GLS`, AIPS spectral forms, `VELREF`, `EPOCH`, `RADECSYS`, and reordered celestial `CROTA`; do not silently repair malformed headers.
- [ ] **P3** Add an opt-in equivalent of WCSLIB's `cylfix` for cylindrical all-sky maps whose pixels extend beyond native longitude ±180°. Astropy does not apply it by default.
- [ ] **P3** Add per-axis invalid-coordinate status or masks together with batch evaluation, so independent valid axes survive failures in another coordinate group. Requires the batch result and error/mask API to be designed first.
- [x] **P2** Expand differential coverage at projection domains, longitude wrapping, cube/facet boundaries, strong projection parameters (including ZPN's first monotonic radial interval), and strong SIP/TPD inverse cases.
- [x] **P3** Expand real-container tests for paired/transposed lookup arrays, alternate-key and `EXTVER` isolation, duplicate matching extensions, and coupled multidimensional `-TAB` tables.

## Other declarative and text formats

- [ ] **P1** Define a stable, human-editable GWCS.jl schema that can be represented in YAML, JSON, or TOML without Python class names.
- [ ] **P1** Support loading the same declarative schema from a string, file, IO object, or already-parsed dictionary.
- [ ] **P1** Add serialization to the declarative schema for reproducible configuration and tests.
- [ ] **P1** Include schema identifiers and versions so files can be migrated safely.
- [ ] **P1** Produce path-aware errors that collect multiple validation failures where possible.
- [ ] **P2** Add a compact line-oriented or table-based format for simple affine and polynomial WCS calibrations.
- [ ] **P2** Add import adapters for observatory calibration tables when a stable external format is identified.

## Performance, compilation, and extensibility

- [ ] **P1** Preserve zero-allocation scalar evaluation as new transforms and frame features are added.
- [ ] **P1** Keep unit conversion, high-level coordinate construction, and metadata bookkeeping outside the compiled numeric kernel.
- [x] **P1** Add benchmark coverage for selectors.
- [x] **P1** Add benchmark coverage for the zenithal projection family.
- [x] **P1** Add benchmark coverage for the remaining projection families.
- [x] **P1** Add benchmark coverage for FITS lookup tables, tabular interpolation, nonlinear spectral algorithms, and grisms.
- [ ] **P1** Add benchmark coverage for mixed-axis frames and batch inversion.
- [x] **P1** Add construction benchmarks for mutation/recompilation.
- [ ] **P1** Add construction benchmarks for large metadata-defined pipelines.
- [x] **P1** Define a stable interface for transforms to provide optimized compilation, inversion, Jacobian, and separability methods.
- [ ] **P2** Add generated or tuple-recursive kernels for large static pipelines when normal compiler inference becomes insufficient.
- [ ] **P2** Investigate CPU-vectorized and threaded batch evaluation without changing scalar semantics.
- [ ] **P2** Investigate GPU-compatible compiled transforms for array workloads where all leaf transforms support the target backend.

## Validation and compatibility coverage

- [ ] **P0** Build a corpus of Python GWCS ASDF examples covering every supported frame and transform tag.
- [x] **P0** Cross-check the zenithal projection family against Astropy/wcslib and FITSWCS.jl.
- [x] **P0** Cross-check each remaining newly supported projection against Astropy/wcslib and FITSWCS.jl where an equivalent exists.
- [ ] **P1** Add round-trip ASDF serialization tests between Julia and Python.
- [ ] **P1** Add FITS import/export round-trip tests with measured approximation errors.
- [x] **P1** Add selector and compound-bounding-box fixtures from IFU-style pipelines.
- [x] **P1** Add mixed celestial/spectral/time frame fixtures and axis-order permutations.
- [ ] **P1** Add randomized forward/inverse property tests for transforms with analytical inverses.
- [ ] **P1** Add numerical-inverse stress tests near projection boundaries, longitude wraps, poles, and highly distorted regions.
- [x] **P1** Add compatibility tests for the oldest supported Julia release and current stable Julia.

## Python GWCS features or design choices we should not duplicate

- **Do not reproduce the Astropy `Model` class hierarchy.** Julia's `AbstractTransform` interface, multiple dispatch, and explicit compilation hooks are simpler and better suited to type-stable scalar kernels. We should implement useful transform behavior, not Astropy's parameter descriptors, model sets, fitter flags, or Python inheritance protocol.
- **Do not implement every Astropy model merely because Python GWCS can compose it.** Source-profile and general modeling functions such as Gaussian, Sersic, Moffat, blackbody, and power-law models are not normally coordinate transformations. External packages can wrap them through the transform-registration API if a real WCS use case appears.
- **Do not duplicate APE 14's Python-specific inheritance and object-description protocol verbatim.** Accessors for dimensions, names, units, physical types, bounds, and correlations are useful, but `world_axis_object_classes`, Python class import strings, and `serialized_classes` exist primarily for Astropy interoperability.
- **Do not recreate `SkyCoord`, `SpectralCoord`, `Time`, or `StokesCoord` inside GWCS.jl.** Optional extensions should interoperate with established Julia packages. The numeric WCS core should remain independent of any one high-level coordinate representation.
- **Do not thread dynamically interpreted units through every scalar transform.** Unitful inputs should be validated and converted at pipeline boundaries, while the compiled inner kernel remains unit-free and allocation-free.
- **Do not expose an unchecked mutable pipeline vector.** Python GWCS acknowledges that direct access can bypass validation. GWCS.jl should use a checked builder or persistent edit operations and produce an immutable compiled WCS.
- **Do not port `EmptyFrame` or legacy string-frame compatibility.** New Julia pipelines should require explicit frame dimensions and metadata, avoiding deprecated behavior that guesses frame structure from transforms.
- **Do not copy Python GWCS's deprecated aliases, indexing shims, warning classes, or historical serialization quirks unless they are required to read real files.**
- **Do not copy Python GWCS's exact fixed-point numerical-inverse algorithm or SciPy fallback.** The ForwardDiff-based damped Newton solver is more idiomatic in Julia. We should add the useful batch diagnostics and initialization strategies while retaining a Julia-native solver architecture.
- **Do not copy the legacy polygon scan-line implementation.** Region support should use a robust Julia geometry or rasterization implementation with clearly defined boundary semantics.
- **Do not reimplement FITS or ASDF container libraries inside GWCS.jl.** FITSIO.jl and FITSFiles.jl should provide FITS container access, and ASDF container handling should delegate to ASDF.jl. GWCS.jl should own FITS WCS interpretation and conversion into compiled pipelines.
- **Do not mirror Astropy's exact NumPy broadcasting and scalar-squeezing rules.** Julia APIs should use normal broadcasting, arrays, tuples, and static vectors while remaining explicit about coordinate-axis layout.
- **Do not make observatory-specific optics and selector machinery mandatory core dependencies.** Grating, Sellmeier, IFU-region, and similar features should be lightweight internal modules or optional extensions unless they become broadly used.
