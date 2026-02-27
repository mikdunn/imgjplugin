# imgjplugin

FFT orientation analysis plugin for Fiji / ImageJ (ImageJ1 "legacy" plugin).

## What it does

The plugin `fftanalysis.imagej.FIBA_Orientation` implements a MATLAB-modeled pipeline (ported from `fiba.m`):

- Contrast stretch (imadjust-like)
- Tukey windowing
- 2D FFT + shifted magnitude display
- Orientation signal extraction $SOL(\theta)$ (0–179°)
- Statistically significant peak detection + weighted-median peak angle
- Band-limited inverse FFT reconstruction + overlay

It produces the same *types* of outputs as the MATLAB script:

- `*_rec.jpg`: composite panel (original, windowed, FFT display, reconstruction mask/result, overlay)
- `*_dat.jpg`: SOL plot with mean and mean+std threshold
- ImageJ **Results** table metrics (`pAng_deg`, `spWid_deg`, `bandStrength`, `pWidth_deg`, `warnPk`, `ang1_deg`, `ang2_deg`)

## Install (Fiji)

1. Build the jar (see below), then copy the shaded jar `*-shaded.jar` into your Fiji `plugins/` folder.
2. Restart Fiji.
3. Run via: `Plugins > FFT > FIBA Orientation (MATLAB)`

## Install (ImageJ2)

ImageJ2 can run ImageJ1 plugins via **legacy** support.

1. Ensure ImageJ2 legacy is enabled/available.
2. Copy the shaded jar `*-shaded.jar` into the appropriate `plugins/` folder for your ImageJ2 installation.
3. Restart ImageJ2.
4. Look for: `Plugins > FFT > FIBA Orientation (MATLAB)`

## Build / smoke test

- Maven project lives in `FFT/`.
- Unit tests live in `FFT/src/test/java` and validate the core math without launching the UI.

This repo also contains a GitHub Actions workflow which builds and uploads the plugin jar as an artifact on each push/PR.

### Note on VS Code + Maven

If you opened this repo using the VS Code **GitHub Repositories** virtual filesystem (`vscode-vfs://...`), your local terminal is *not* running inside a real folder on disk. In that case:

- `cd FFT` will fail (because there is no `FFT/` directory on your C: drive)
- `mvn` will fail unless Maven is installed and on your `PATH`

Two ways to build:

1. Use the GitHub Actions artifact (no local Maven needed):
	- Commit + push your changes (a push to GitHub triggers the workflow).
	- In GitHub, go to **Actions** → workflow **build** → open the latest run.
	- Under **Artifacts**, download `imgjplugin-fft`.
	- Unzip it and use the `FFT/target/*-shaded.jar` inside.
2. Clone the repo locally (so it exists on disk), install Maven, then run `mvn test package` in `FFT/`.

### Running in Fiji / getting output images

1. Copy the `*-shaded.jar` into your Fiji `plugins/` folder.
2. Restart Fiji.
3. Run the plugin from `Plugins > FFT > FIBA Orientation (MATLAB)`.
4. The plugin writes outputs alongside the input image (filenames like `*_rec.jpg` and `*_dat.jpg`).

### One-click macro runner

If you don’t want to click through menus each time, use the included macro:

- Run `macros/run_fiba_orientation_outputs.ijm` via `Plugins > Macros > Run...`
- It runs `FIBA Orientation (MATLAB)` on the currently active image and saves the standard outputs.

## Notebook environment preflight (timeout prevention)

`notebooks/fiba_tile_montage_pipeline.ipynb` now includes a **preflight cell** that checks:

- Required commands on `PATH`: `python`, `java`, `quarto`, `jupyter`, `mvn`
- Required Python packages: `numpy`, `matplotlib`, `Pillow`, `tifffile`

Use the setup script before running the notebook pipeline cells:

- `notebooks/setup_notebook_env.ps1`
- `notebooks/requirements-fiba-notebook.txt`

The script upgrades/install Python dependencies and adds the Python user `Scripts` folder to your user `PATH` (so `jupyter` is discoverable).
If `mvn` is still missing after setup, install Maven and add its `bin` directory to `PATH`.
