// Run the FFT orientation plugin on the currently active image.
// Produces outputs (`*_rec.jpg`, `*_dat.jpg`) next to the input image (or use outputdir=...).
//
// Usage (in Fiji/ImageJ): Plugins > Macros > Run... and select this file.

requires("1.52p");

// Ensure there is an active image.
if (nImages==0) {
    exit("No image is open. Open an image first, then run this macro.");
}

// Run the plugin. Options are parsed as whitespace-separated key=value pairs.
// Artifact suppression is intentionally NOT enabled here; we will deal with artifacts ONLY with alpha values.
run("FIBA Orientation (MATLAB)",
    "save=true showcomposite=true showplot=true " +
    "");
