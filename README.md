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

1. Build the jar (see below), then copy the shaded jar `*-all.jar` into your Fiji `plugins/` folder.
2. Restart Fiji.
3. Run via: `Plugins > FFT > FIBA Orientation (MATLAB)`

## Install (ImageJ2)

ImageJ2 can run ImageJ1 plugins via **legacy** support.

1. Ensure ImageJ2 legacy is enabled/available.
2. Copy the shaded jar `*-all.jar` into the appropriate `plugins/` folder for your ImageJ2 installation.
3. Restart ImageJ2.
4. Look for: `Plugins > FFT > FIBA Orientation (MATLAB)`

## Build / smoke test

- Maven project lives in `FFT/`.
- Unit tests live in `FFT/src/test/java` and validate the core math without launching the UI.

This repo also contains a GitHub Actions workflow which builds and uploads the plugin jar as an artifact on each push/PR.