# Analytic Magnetic Helical Field — Cross-Platform GUI

This package wraps the analytic C++ magnetic-field solver and Python plotting pipeline in a desktop GUI. The source version now supports Linux, macOS, and native Windows through the same CMake build path, and it can also be frozen into standalone applications that do not require the end user to install Python, Boost, CMake, or a C++ compiler.

## Application workflow

The GUI has three tabs:

1. **Setup**
   - Edit the physical and numerical solver parameters.
   - See a live 3D preview of all `N` helical strands, their centerlines, the outer cylindrical wall, and the **true intersections of each strand surface with its helix-normal plane**. The reference strand's plane is also shown as a translucent guide. The strands are evenly spaced by azimuthal phases `2π j/N`.
   - The geometry updates as parameters are edited.

2. **Run**
   - In source mode, configures and builds `main.cpp` with CMake when needed.
   - In a packaged standalone app, uses the native solver bundled inside the application; the user's machine does not need a compiler.
   - Writes the GUI settings to a run-specific `config.txt`.
   - Runs the C++ calculation and displays solver progress and console output.
   - Runs plotting in a separate worker process, then loads the results automatically.

3. **Results**
   - **Plots** sub-tab: previews every generated PNG, saves the selected plot, or exports all plots.
   - **CSV data** sub-tab: previews the user-facing tabular CSV outputs in a scrollable table and can save the selected CSV or all CSVs to any directory. The preview is intentionally limited to the first 250 data rows so large field files do not freeze the GUI.
   - `Save all CSVs…` includes the raw compact solver `field.csv`; there is also a dedicated `Save raw solver field.csv…` button.
   - Shows numerical diagnostics from `summary.json`.
   - The full compact solver output is not opened as a spreadsheet preview because its five logical rows can each contain enormous numbers of columns; the Results table uses the row-oriented processed CSVs described below.

Runs are stored under `~/AnalyticMagneticField/runs/` by default, so a standalone application's temporary installation directory is never used for user data.

---

# Requirements for running from source

The Python requirements are in `requirements.txt`:

- NumPy
- SciPy
- Matplotlib
- Pillow
- Tkinter from your Python installation/OS

The native solver additionally requires:

- a C++20 compiler
- CMake 3.16 or newer
- Boost headers

Boost is header-only for this program; no Boost libraries need to be linked.

The old Makefile is retained for convenience on Unix-like systems, but the GUI itself now uses **CMake** on every operating system.

---

# Linux

## Ubuntu / Debian

Install system dependencies:

```bash
sudo apt update
sudo apt install build-essential cmake libboost-all-dev python3-tk python3-venv
```

Create the Python environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Launch:

```bash
./launch.sh
```

or:

```bash
python3 app.py
```

The GUI will configure/build the C++ solver automatically with CMake when you press **Run simulation**.

You can also build it manually:

```bash
python3 build_solver.py
```

---

# macOS

Install Xcode command-line tools:

```bash
xcode-select --install
```

Install CMake and Boost. With Homebrew:

```bash
brew install cmake boost
```

Use a Python installation with Tk support. Then create a virtual environment:

```bash
python3 -m venv .venv
source .venv/bin/activate
pip install -r requirements.txt
```

Launch with either:

```bash
./launch.command
```

or:

```bash
python3 app.py
```

If CMake does not find a Homebrew Boost installation automatically, enter the Homebrew prefix in **Boost include root (optional)** in the Setup tab. Typical examples are:

```text
/opt/homebrew
```

on Apple Silicon, or:

```text
/usr/local
```

on older Intel/Homebrew installations.

Unlike the earlier version of this project, GNU GCC/OpenMP is no longer required: the solver does not currently use OpenMP pragmas, so Apple Clang works normally.

---

# Windows 10 / 11 — native

The application no longer requires WSL. It can build and run natively on Windows.

Install:

1. **Python 3** from python.org. The standard Windows installer includes Tkinter.
2. **CMake**, with the option to add CMake to `PATH`.
3. **Visual Studio 2022 Build Tools** with the **Desktop development with C++** workload.
4. **Boost headers**. You may simply download and extract a Boost release; building Boost itself is unnecessary for this project.

For example, if Boost is unpacked as:

```text
C:\Libraries\boost_1_90_0\boost\...
```

then use:

```text
C:\Libraries\boost_1_90_0
```

as **Boost include root (optional)** in the GUI.

Create an environment from PowerShell or Command Prompt:

```bat
py -3 -m venv .venv
.venv\Scripts\activate
pip install -r requirements.txt
```

Launch by double-clicking `launch.bat`, or run:

```bat
launch.bat
```

You may also use PowerShell:

```powershell
.\launch.ps1
```

or directly:

```bat
py -3 app.py
```

When Run is pressed, CMake will use the installed Visual Studio generator and build `field_solver.exe` automatically.

The CMake configuration uses the static MSVC runtime for the native solver, reducing deployment dependencies for packaged Windows builds.

---

# Building the solver manually on any OS

The helper script is cross-platform:

```bash
python build_solver.py
```

For a non-system Boost installation:

```bash
python build_solver.py --boost-root /path/to/boost
```

Windows example:

```bat
py -3 build_solver.py --boost-root C:\Libraries\boost_1_90_0
```

To discard the previous CMake build first:

```bash
python build_solver.py --clean
```

Equivalent raw CMake commands are:

```bash
cmake -S . -B build -DCMAKE_BUILD_TYPE=Release
cmake --build build --config Release --parallel
```

The executable is normally one of:

```text
Linux/macOS: build/field_solver
Windows:     build/Release/field_solver.exe
```

The application searches both single-configuration and multi-configuration CMake layouts automatically.

---

# Creating standalone applications

Standalone builds use **PyInstaller**. A packaged user does **not** need Python, CMake, Boost, or a compiler; the package contains the Python runtime/dependencies and the already-compiled native C++ solver.

## Important: build separately on each OS

PyInstaller is not a cross-compiler. Build the Windows application on Windows, the macOS application on macOS, and the Linux application on Linux.

Install packaging dependencies in your virtual environment:

```bash
pip install -r requirements-packaging.txt
```

Then run:

```bash
python package_app.py --clean
```

This performs two steps automatically:

1. builds the C++ solver in Release mode with CMake;
2. invokes PyInstaller and bundles that native solver inside the GUI application.

The default is **one-directory** mode because it starts faster and is easier to debug for a NumPy/SciPy/Matplotlib application.

To request a single-file build instead:

```bash
python package_app.py --mode onefile --clean
```

If Boost is not installed system-wide:

```bash
python package_app.py --boost-root /path/to/boost --clean
```

or on Windows:

```bat
py -3 package_app.py --boost-root C:\Libraries\boost_1_90_0 --clean
```

## Windows standalone result

Recommended build command:

```bat
py -3 package_app.py --clean
```

The output is under:

```text
dist\AnalyticMagneticField\
```

The user launches:

```text
AnalyticMagneticField.exe
```

For a single `.exe` instead:

```bat
py -3 package_app.py --mode onefile --clean
```

which produces approximately:

```text
dist\AnalyticMagneticField.exe
```

A one-directory build is recommended first; once that has been tested on a clean Windows machine, a one-file build is reasonable for easier distribution.

## macOS standalone result

Build on the Mac architecture you intend to support:

```bash
python package_app.py --clean
```

PyInstaller's windowed macOS build produces an application bundle under `dist`, normally:

```text
dist/AnalyticMagneticField.app
```

You can launch it in Finder like any other application.

For distribution to other users outside your own Mac, Apple Gatekeeper generally requires **code signing**, and broad public distribution normally also requires **notarization**. Those are publishing steps after the application bundle has been built and tested.

Apple Silicon and Intel builds are architecture-specific unless you deliberately build with a universal2 Python/toolchain. Build/test on each architecture you intend to distribute, or use a universal2 environment.

## Linux standalone result

Build:

```bash
python package_app.py --clean
```

The bundle is under:

```text
dist/AnalyticMagneticField/
```

and the executable inside that directory is `AnalyticMagneticField`.

Linux binary compatibility depends on the system C library. For the broadest compatibility, build the release on the **oldest Linux distribution you intend to support**, then test it on newer distributions.

For a single binary:

```bash
python package_app.py --mode onefile --clean
```

## Why the native solver is bundled

In source mode, the GUI invokes CMake because developers may change `main.cpp`.

In standalone mode, recompiling on an end user's computer would defeat the purpose of a standalone application. `package_app.py` therefore compiles `field_solver` during packaging and adds it under the application's internal `bin/` directory. At runtime the GUI detects the frozen package and launches the bundled solver directly.

The plotting stage is also packaging-aware. A frozen GUI cannot safely assume `sys.executable` is a normal Python interpreter, so the application re-launches itself in a private `--plot-worker` mode. This keeps Matplotlib plotting isolated from the live Tk GUI and works in both source and frozen modes.

---

# Parameters exposed in the GUI

- `I` — total current [A]
- `a` — helix center radius [m]
- `b` — strand radius [m]
- `wavevector` — helical wavevector `k` [m^-1]
- `outer_r` — outer cylindrical wall radius [m]
- `delta` — radial integration/sampling step [m]
- `N` — **number of helical strands**. The GUI and generated geometry plots draw all `N` strands, evenly spaced in azimuth by `2π/N`. Because the modeled strand surfaces are screw-rotated circular cross-sections, neighboring strand centers are separated at every fixed `z` by the exact chord `d = 2 a sin(π/N)`. The GUI warns before a run if the strand surfaces touch or overlap (`d <= 2b`).
- `index_cap` — maximum Fourier index
- `amb_thetas` — base angular sampling count
- `contour_points` — resolution of the plane/helix-surface intersection used in post-processing
- `plot_dpi` — saved-plot resolution
- `xy_grid_resolution` — number of Cartesian samples per axis for the smooth general `x-y` field plots and `field_xy_plane.csv` export
- `preview_turns` — number of helix turns shown in the live preview only

The original constraints are validated before a run, including `0 < b < a`, `outer_r > a + b`, and a nonzero wavevector.

---

# Generated run files

A normal run directory contains:

- `config.txt` — exact solver configuration
- `gui_settings.json` — GUI/post-processing settings
- `field.csv` — full native solver output in the compact five-line polar format. This is the authoritative raw field data used by post-processing.
- `field_xy_plane.csv` — smooth, row-oriented Cartesian `z=0` field grid used for the general `x-y` plots. Columns include grid indices, `x,y,r,theta`, `B_x,B_y,B_z`, in-plane magnitude, and total magnitude. This is the convenient general-field CSV shown in the GUI.
- `b_normal_plane_contour.csv` — magnetic field sampled on the reference strand's true helix-normal plane/surface intersection. This is also the data source for the contour-field plots.
- `normal_plane_contours_geometry.csv` — geometry of the corresponding new contour on every one of the `N` strands.
- `field_xy_magnitude.png` — total magnetic-field magnitude across the general `x-y` plane.
- `field_xy_bz.png` — `B_z` across the general `x-y` plane.
- `field_xy_inplane.png` — in-plane field magnitude with `B_x,B_y` direction arrows.
- `b_normal_plane_contour.png` — total field magnitude, normal component, and in-plane magnitude around the new helix-normal contour, plotted against contour arc length.
- `b_normal_plane_contour_components.png` — Cartesian `B_x`, `B_y`, and `B_z` components around the same contour.
- `normal_plane_contours_on_helix_3d.png` — geometry-only validation view showing all `N` new contours on the full helical geometry.
- `summary.json`
- `plots_manifest.json`
- `csv_manifest.json` — user-facing tabular CSV files exposed in the Results/CSV-data preview.

---

# Plane/surface contour

For the right-hand strand centerline

```text
c(z) = (a cos(kz), a sin(kz), z),
```

the local tangent at `z=0` is

```text
(0, a k, 1).
```

For strand 0, the plotting pipeline therefore intersects the exact helical strand surface with the fixed plane

```text
a k y + z = 0.
```

The closed branch is followed numerically around the surface and checked for closure and plane-equation residual. For strand `j`, the strand and its local normal plane are both obtained by the azimuthal rotation `2π j/N`, so its new contour is the same rigid rotation of the reference contour. The GUI preview and the generated 3D plot use this same geometry routine.

The old `z=0` circular strand-intersection contours and their dedicated field-around-the-contour plots are no longer generated. The application does, however, generate general magnetic-field maps over the entire `x-y` plane.

---

# Ripple fix / interpolation

The original contour lookup selected the nearest discrete solver sample, producing visible small jumps as the selected index changed. The application instead interpolates on the solver's native polar rings: angular interpolation is periodic, and neighboring radial rings are linearly interpolated.

The C++ output writer uses high-precision floating-point output so coordinate and field values are not unnecessarily quantized before post-processing.

---

# Distribution checklist

Before giving a standalone package to other users:

1. Build it on the same operating system family you are targeting.
2. Prefer `--mode onedir` for the first release and test the resulting folder on a clean machine.
3. Run a small calculation and verify the three general `x-y` field plots plus the 3D contour-geometry plot are generated.
4. Test both Results sub-tabs: plot preview/export and CSV preview/export.
5. Confirm `field.csv`, `field_xy_plane.csv`, and `b_normal_plane_contour.csv` are present in the run directory.
6. On macOS, sign/notarize the `.app` before broad external distribution.
7. On Linux, build on the oldest supported distribution for better `glibc` compatibility.
8. On Windows, test on a machine without your development environment to verify the package is truly self-contained.

